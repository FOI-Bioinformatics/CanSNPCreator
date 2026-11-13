"""Defects found by reviewing the code, each with the input that showed it."""

import gzip
import logging
from pathlib import Path

import pytest
from conftest import cli_run

from cansnpscorer.helpers.data_utils import codon_complement, search
from cansnpscorer.helpers.file_utils import (
    read_fasta,
    read_genome_map,
    read_intervals,
)
from cansnpscorer.helpers.names import resolve_clusters
from cansnpscorer.hierarchy import Hierarchy, HierarchyError
from cansnpscorer.model import Reference
from cansnpscorer.scheme import Scheme, SchemeSnp, hierarchy_from_paths, schemes_in

# -- placing genomes and auditing ----------------------------------------------


def scheme_of(paths, snps):
    hierarchy, _conflicts, _problems = hierarchy_from_paths(paths)
    rows = [SchemeSnp(node=node, strain="s", reference="r", genome=genome, position=position,
                      derived=derived, ancestral=ancestral)
            for node, genome, position, derived, ancestral in snps]
    return Scheme(organism="t", paths=paths, snps=rows, hierarchy=hierarchy)


def test_place_refuses_a_scheme_on_two_genomes_without_genome():
    # It took whichever genome's SNP came first for each node.
    from cansnpscorer.tree_tools import place

    scheme = scheme_of([["R"], ["R", "N"]], [("N", "G2", 3, "T", "G"), ("N", "G1", 3, "T", "G")])
    reference = Reference.from_records([("c1", "ACGTACGT")])
    alignment = {"r": "ACGTACGT", "a": "ACTTACGT", "b": "ACTTACGT"}
    tree = Hierarchy.from_newick("(r,(a,b));", support_labels=True)
    with pytest.raises(ValueError, match="--genome must say"):
        place(scheme, tree, alignment, reference)
    with pytest.raises(ValueError, match="no SNP on genome 'G9'"):
        place(scheme, tree, alignment, reference, genome="G9")
    assert place(scheme, tree, alignment, reference, genome="G1").clusters


def test_a_parent_holding_its_childs_genomes_is_offered_the_childs_site():
    # P has the single child C; both hold a and b. P's own SNP is carried by
    # a alone, and another site separates exactly a and b. That site fits P as
    # well as C, so P is offered a replacement, not called "not a clade".
    from cansnpscorer.alignment_audit import check_against_alignment

    scheme = scheme_of([["R"], ["R", "P"], ["R", "P", "C"], ["R", "X"]],
                       [("P", "G", 2, "T", "C")])
    reference = Reference.from_records([("c1", "ACGTACGT")])
    alignment = {"r": "ACGTACGT", "a": "ATGTACGA", "b": "ACGTACGA", "x": "ACGTACGT"}
    findings = check_against_alignment(scheme, alignment, reference, {"C": ["a", "b"],
                                                                      "X": ["x"]})
    categories = {finding.category for finding in findings}
    assert "unsupported-node" not in categories
    assert "replacement" in categories


def test_audit_from_the_genomes_file_keeps_the_genome_option(tmp_path, caplog):
    # A loop over the genomes file reused the name genome, which then held the
    # last accession and made every SNP of the scheme be skipped.
    from cansnpscorer.cansnp_audit import main as audit_main

    directory = tmp_path / "scheme"
    directory.mkdir()
    (directory / "t_tree.txt").write_text("R\nR\tN\n")
    (directory / "t_snp.txt").write_text(
        "snp_id\tstrain\treference\tgenome\tposition\tderived_base\tancestral_base\n"
        "R\ts\tr\tREF\t1\tA\tG\nN\ts\tr\tREF\t3\tT\tG\n")
    (directory / "t_references.txt").write_text(
        "genome\tstrain\tgenbank_id\trefseq_id\tassembly_name\n"
        "REF\ts\tGCA_000000001.1\tGCF_000000001.1\tASM1v1\n")
    (directory / "t_genomes.txt").write_text(
        "genome\tnode\trole\nGCA_000000002.1\tN\trepresentative\n"
        "GCA_000000003.1\tN\trepresentative\n")
    reference = tmp_path / "ref.fasta"
    reference.write_text(">c1\nACGTACGT\n")
    alignment = tmp_path / "aln.fasta"
    # GCA_000000003.1 lacks N's derived base T.
    alignment.write_text(">x_GCA_000000001.1\nACGTACGT\n>a_GCA_000000002.1\nACTTACGT\n"
                         ">b_GCA_000000003.1\nACGTACGT\n")
    report = tmp_path / "findings.tsv"
    audit_main(str(directory), reference=str(reference), alignment=str(alignment),
               output=str(report))
    assert "not-canonical" in report.read_text()


def test_audit_refuses_a_genome_the_scheme_does_not_name(written, dataset, caplog):
    from cansnpscorer.cansnp_audit import main as audit_main

    outdir, _scores = written
    assert audit_main(outdir, reference=dataset.reference, genome="ZZ") == 1
    assert "no SNP on genome 'ZZ'" in caplog.text


def test_place_leaves_single_leaves_out_as_clades():
    # On an unresolved tree the genomes outside one leaf stood in for a clade,
    # and d, which lacks the derived base, was placed at N1.
    from cansnpscorer.tree_tools import place

    scheme = scheme_of([["R"], ["R", "N1"]], [("N1", "G", 3, "T", "G")])
    reference = Reference.from_records([("c1", "ACGTACGT")])
    alignment = {"r": "ACGTACGT", "a": "ACTTACGT", "b": "ACTTACGT", "c": "ACGTACGT",
                 "d": "ACGTACGT"}
    # The tree holds the four genomes but not the reference, as the reviewer's
    # star of four did: outside c is then {a, b, d}, which matched N1 at 0.67.
    placement = place(scheme, Hierarchy.from_newick("(a,b,c,d);", support_labels=True),
                      alignment, reference)
    assert "d" not in placement.clusters.get("N1", [])


def test_node_table_supports_in_percent_are_read_as_fractions(tmp_path):
    from cansnpscorer.tree_tools import collapse, read_genome_tree

    edges = tmp_path / "tree2tax.tsv"
    edges.write_text("child\tparent\nn1\troot\nn2\troot\na\tn1\nb\tn1\nc\tn2\nd\tn2\n")
    nodes = tmp_path / "nodes.tsv"
    nodes.write_text("node\tparent\tn_leaves\tdepth\tbranch_length\tsupport\tred\tlineage\n"
                     "n1\troot\t2\t1\t1\t95\tNA\t\nn2\troot\t2\t1\t1\t30\tNA\t\n")
    tree = read_genome_tree(str(edges), str(nodes))
    assert tree.support == {"n1": 0.95, "n2": 0.30}
    assert set(collapse(tree, min_support=0.9).removed) == {"n2"}


def test_a_member_under_two_leaves_keeps_the_first(tmp_path, caplog):
    from cansnpscorer.tree_tools import read_members

    genome_map = tmp_path / "genomes_map.tsv"
    genome_map.write_text("m1\ta\nm1\tc\nd\ta\n")
    members = read_members(str(genome_map), Hierarchy.from_newick("((a,b),(c,d));"))
    assert members == {"a": ["m1"]}, "m1 kept under a; d is a leaf itself"
    assert "under two leaves" in caplog.text


def test_place_reports_an_empty_genome_map(written, dataset, tmp_path):
    # A tree whose leaves are none of the alignment's genomes places nothing;
    # that used to exit 0 with an empty map.
    outdir, _scores = written
    tree = tmp_path / "tree.nwk"
    tree.write_text("((x1,x2),(x3,x4));")
    assert cli_run("place", "--scheme-dir", outdir, "--tree", str(tree),
                   "--alignment", dataset.alignment, "--reference", dataset.reference,
                   "-o", str(tmp_path / "placed.tsv")) == 2


# -- trees ---------------------------------------------------------------------


def test_a_generated_name_avoids_a_real_label():
    # The unnamed node took the name node1, which a real node already had.
    tree = Hierarchy.from_newick("((a,b),(c,d)node1);")
    assert sorted(tree.children["node1"]) == ["c", "d"]
    assert len(tree.leaves) == 4


def test_a_name_used_twice_is_refused():
    with pytest.raises(HierarchyError, match="names 'a' twice"):
        Hierarchy.from_newick("(a,a,b);")


def test_a_quoted_label_may_hold_brackets_and_colons():
    tree = Hierarchy.from_newick("('(x)',':1')r;")
    assert set(tree.leaves) == {"(x)", ":1"}


# -- scoring -------------------------------------------------------------------


def test_search_finds_the_nearest_interval():
    # The binary search stopped at the first window holding the position.
    assert search(12, [(0, 10), (30, 40), (100, 110)], 25) == (pytest.approx(0.08), None)
    assert search(35, [(0, 10), (30, 40), (100, 110)], 25) == (0.0, (30, 40))
    assert search(70, [(0, 10), (30, 40), (100, 110)], 25) == (1.0, None)


def test_masks_from_several_names_are_joined_on_their_contig():
    from cansnpscorer.cansnp_evaluator import in_regions, regions_over_contigs

    placed = regions_over_contigs({"x": [(0, 100)], "y": [(10, 20)]},
                                  Reference.from_records([("c1", "A" * 200)]))
    assert placed == {"c1": [(0, 100)]}
    assert in_regions(50, placed["c1"])


def test_codon_complement_reads_n():
    assert codon_complement("ANG") == "TNC"


# -- files ---------------------------------------------------------------------


def test_bed_contigs_named_like_track_lines_are_kept(tmp_path):
    bed = tmp_path / "r.bed"
    bed.write_text("track name=x\ntrackA\t1\t5\nbrowserX\t2\t4\nchr1 0 3\n")
    assert read_intervals(str(bed)) == {"trackA": [(1, 4)], "browserX": [(2, 3)],
                                        "chr1": [(0, 2)]}


def test_gff3_after_a_comment_or_byte_order_mark_is_read_as_gff3(tmp_path):
    gff = tmp_path / "r.gff"
    gff.write_bytes(b"\xef\xbb\xbf# made by hand\n##gff-version 3\n"
                    b"c1\tx\tregion\t10\t20\t.\t+\t.\tID=a\n")
    assert read_intervals(str(gff)) == {"c1": [(9, 19)]}


def test_fasta_comments_spaces_and_bad_characters(tmp_path):
    good = tmp_path / "a.fasta"
    good.write_text(">a description\nAC GT\n;a comment\nNN\n")
    assert read_fasta(str(good)) == [("a", "ACGTNN")]
    for text, message in [(">a\nAC1T\n", "not IUPAC"), (">a\n>b\nAC\n", "hold no sequence"),
                          ("> a\nAC\n", "header without an ID")]:
        bad = tmp_path / "b.fasta"
        bad.write_text(text)
        with pytest.raises(ValueError, match=message):
            read_fasta(str(bad))


def test_a_truncated_gzip_file_is_an_error_not_a_traceback(dataset, tmp_path, caplog):
    whole = gzip.compress(Path(dataset.reference).read_bytes())
    truncated = tmp_path / "reference.fasta.gz"
    truncated.write_bytes(whole[: len(whole) // 2])
    with caplog.at_level(logging.ERROR):
        status = cli_run("trf", "--reference", str(truncated), "-o", str(tmp_path / "r.bed"))
    assert status == 1
    status = cli_run("extract", "--alignment", dataset.alignment,
                     "--reference", str(truncated), "--hierarchy", dataset.hierarchy,
                     "--genome-map", dataset.genome_map, "-o", str(tmp_path / "s.vcf"))
    assert status == 1


def test_a_genome_map_header_is_not_a_node(tmp_path):
    genome_map = tmp_path / "map.tsv"
    genome_map.write_text("genome\tnode\ng1\tA\n")
    assert read_genome_map(str(genome_map)) == {"A": ["g1"]}


def test_two_names_of_one_record_at_two_nodes_are_refused():
    resolved, _unmatched, ambiguous = resolve_clusters({"n1": ["A_1"], "n2": ["A.1"]}, ["A-1"])
    assert sorted(ambiguous) == ["A.1", "A_1"]
    assert resolved == {"n1": [], "n2": []}


def test_vcf_contig_ids_and_spaces_follow_vcf_4_3(tmp_path):
    from cansnpscorer import vcf

    with pytest.raises(vcf.VcfError, match="does not allow as IDs"):
        vcf.header_lines(Reference.from_records([("c,1", "ACGT")]), "r.fasta", "t", [])
    record = vcf.snp_record(Reference.from_records([("c1", "ACGT")]), "c1", 1, ("A", "G"),
                            node="Clade 1", ancestral="A")
    assert "NODE=Clade%201" in record.line()


def test_a_scheme_directory_with_glob_characters_is_found(tmp_path):
    directory = tmp_path / "d[1]"
    directory.mkdir()
    (directory / "x_tree.txt").write_text("R\n")
    assert schemes_in(str(directory)) == ["x"]


@pytest.mark.parametrize("field,value,message", [
    ("genbank_id", "GCA_000007125.1\n", "well formed GenBank"),
    ("strain", "a\tb", "tab or a line break"),
    ("organism", "a/b", "names the files"),
])
def test_export_metadata_that_would_corrupt_the_files_is_refused(field, value, message):
    from cansnpscorer.cansnp_export import check_metadata

    values = {"genome": "G", "strain": "s", "reference": "r", "genbank_id": "GCA_000007125.1",
                  "refseq_id": "GCF_000007125.1", "assembly_name": "A", "organism": "t"}
    values[field] = value
    assert any(message in problem for problem in check_metadata(**values))


# -- the command line ----------------------------------------------------------


@pytest.mark.parametrize("arguments", [
    ["--gc-goal", "7"], ["--top-n", "0"], ["--tandem-limit", "-1"],
])
def test_out_of_range_options_are_refused_by_the_parser(dataset, tmp_path, arguments):
    with pytest.raises(SystemExit) as stopped:
        cli_run("evaluate", "--snps", "x.vcf", "--reference", dataset.reference,
                "--genes", dataset.gff, "--repeats", dataset.repeats,
                "-o", str(tmp_path / "s.vcf"), *arguments)
    assert stopped.value.code == 2


def test_an_unwritable_output_is_an_error_not_a_traceback(dataset, tmp_path):
    assert cli_run("extract", "--alignment", dataset.alignment,
                   "--reference", dataset.reference, "--hierarchy", dataset.hierarchy,
                   "--genome-map", dataset.genome_map,
                   "-o", str(tmp_path / "no" / "such" / "dir" / "s.vcf")) == 1
