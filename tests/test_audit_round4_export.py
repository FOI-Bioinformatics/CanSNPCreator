"""Defects of export and the scheme reader found in the fourth audit round, each
pinned by the input that showed it."""

from conftest import METADATA, info_of, vcf_lines

from cansnpscorer.cansnp_export import check_metadata, check_scheme, genome_rows
from cansnpscorer.cansnp_export import main as export_main
from cansnpscorer.hierarchy import Hierarchy
from cansnpscorer.model import Reference
from cansnpscorer.scheme import SchemeSnp, read_snp_set

TREE = Hierarchy([("root", None), ("A", "root"), ("A.1", "A"), ("B", "root")])


def test_a_position_given_at_two_nodes_is_refused(written, dataset, tmp_path, caplog):
    _outdir, scores = written
    header, records = vcf_lines(scores)
    passing = [line.split("\t") for line in records if line.split("\t")[6] == "PASS"]
    first = passing[0]
    other = next(fields for fields in passing
                 if info_of(fields)["NODE"] != info_of(first)["NODE"])
    copy = list(first)
    copy[7] = copy[7].replace("NODE=" + info_of(first)["NODE"],
                              "NODE=" + info_of(other)["NODE"], 1)
    joined = tmp_path / "joined.vcf"
    joined.write_text("\n".join(header + records + ["\t".join(copy)]) + "\n")
    outdir = tmp_path / "out"
    assert export_main(str(joined), dataset.hierarchy, dataset.reference, str(outdir),
                       **METADATA) == 2
    assert f"{first[0]}:{first[1]} at " in caplog.text
    assert "more than one passing record" in caplog.text
    assert not outdir.exists()


def test_a_node_name_with_a_tab_or_line_break_is_refused():
    reference = Reference.from_records([("c1", "ACGT")])
    for name in ("A\tx", "A\nx", "A\rx"):
        tree = Hierarchy([("root", None), (name, "root")])
        problems, _warnings = check_scheme(tree, {}, reference)
        assert any(repr(name) in problem for problem in problems), name


def test_the_representative_of_an_assembly_is_kept_whatever_the_order():
    # The GCA_ name of the assembly comes first and is not aligned; the GCF_
    # name is aligned. Before, the first met was kept, as a member.
    clusters = {"A.1": ["GCA_000000021.1", "GCF_000000021.1"],
                "B": ["GCA_000000022.1"], "A": ["GCF_000000022.1"]}
    aligned = ["ref_GCF_000000001.1", "GCF_000000021.1", "x_GCF_000000022.1"]
    rows, warnings = genome_rows(TREE, TREE.nodes, clusters, aligned)
    assert ("GCF_000000021.1", "A.1", "representative") in rows
    assert ("GCF_000000022.1", "A", "representative") in rows
    assert not any(row[0] in ("GCA_000000021.1", "GCA_000000022.1") for row in rows)
    assert any("GCA_000000022.1/GCF_000000022.1 at A and B" in warning for warning in warnings)


def test_an_unplaced_reference_is_listed_as_an_outgroup_with_a_warning():
    rows, warnings = genome_rows(TREE, TREE.nodes, {"A.1": ["GCA_000000002.1"]},
                                 ["ref_GCF_000000001.1", "GCA_000000002.1"])
    assert ("GCF_000000001.1", "", "outgroup") in rows
    assert any("reference genome 'ref_GCF_000000001.1' is placed at no node" in warning
               for warning in warnings)
    # Placed, it is a representative and nothing is said.
    rows, warnings = genome_rows(TREE, TREE.nodes, {"A": ["ref_GCF_000000001.1"]},
                                 ["ref_GCF_000000001.1"])
    assert rows == [("GCF_000000001.1", "A", "representative")]
    assert warnings == []


def test_accession_versions_have_one_or_two_digits_as_cansnper2_data_requires():
    values = {"genome": "G", "strain": "s", "reference": "r",
              "genbank_id": "GCA_000007125.12", "refseq_id": "GCF_000007125.12",
              "assembly_name": "A", "organism": "t"}
    assert check_metadata(**values) == []
    values["genbank_id"] = "GCA_000007125.123"
    assert any("well formed GenBank" in problem for problem in check_metadata(**values))

    rows, warnings = genome_rows(TREE, TREE.nodes,
                                 {"A": ["GCA_000000003.12"], "B": ["GCA_000000004.123"]},
                                 ["ref"])
    assert [row[0] for row in rows] == ["GCA_000000003.12"]
    assert any("GCA_000000004.123" in warning and "one or two digits" in warning
               for warning in warnings)


def test_ranks_of_a_node_that_skip_or_repeat_are_a_problem(tmp_path):
    path = tmp_path / "t_canSNPs.txt"
    path.write_text(
        "node\trank\tgenome\tposition\tderived_base\tancestral_base\tscore\n"
        "A\t1\tG\t10\tA\tC\t9\nA\t3\tG\t20\tA\tC\t8\n"
        "B\t1\tG\t30\tA\tC\t9\nB\t1\tG\t40\tA\tC\t8\n"
        "R\t1\tG\t50\tA\tC\t9\nR\t2\tG\t60\tA\tC\t8\n")
    snps = [SchemeSnp(node=node, strain="", reference="", genome="G", position=position,
                      derived="A", ancestral="C")
            for node, position in (("A", 10), ("B", 30), ("R", 50))]
    problems: list[str] = []
    read_snp_set(str(path), snps, problems)
    assert any("node 'A' the ranks 1, 3" in problem for problem in problems)
    assert any("node 'B' the ranks 1, 1" in problem for problem in problems)
    assert not any("'R'" in problem for problem in problems)
