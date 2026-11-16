"""Genome trees: reading them, collapsing them into a scheme, and placing
genomes on an existing scheme by the tree rather than by typing."""

import logging
from pathlib import Path

import pytest
from conftest import CLUSTERS, PLANTED, cli_run, vcf_records

from cansnpscorer.cansnp_audit import audit_scheme
from cansnpscorer.helpers.data_utils import canSNP_extractor, snp_extractor
from cansnpscorer.helpers.file_utils import (
    read_alignment,
    read_genome_map,
    read_hierarchy,
    reference_parser,
)
from cansnpscorer.helpers.names import resolve_names
from cansnpscorer.hierarchy import Hierarchy, HierarchyError
from cansnpscorer.model import Reference
from cansnpscorer.scheme import Scheme, SchemeSnp, hierarchy_from_paths, read_scheme
from cansnpscorer.tree_tools import (
    collapse,
    place,
    place_by_clade,
    read_genome_tree,
    read_members,
    rename,
    write_collapsed,
)

#: The planted dataset as a genome tree, rooted on the reference genome.
PLANTED_TREE = "(ref:0.5,(((s1:0.01,s1b:0.01)97:0.2,s2:0.3)99:0.4,(s3:0.2,s4:0.2)100:0.4)100:0.1);"


# -- reading Newick ------------------------------------------------------------


def test_numeric_internal_labels_are_support_and_not_names():
    tree = Hierarchy.from_newick("((a:1,b:2)95:0.5,(c,d)80:0.25,e);", support_labels=True)
    assert sorted(tree.support.values()) == [0.8, 0.95], "percentages read as fractions"
    assert "95" not in tree and "80" not in tree
    assert tree.lengths["a"] == 1.0 and len(tree.lengths) == 4


def test_a_schemes_numeric_node_names_stay_names(tmp_path):
    # A hierarchy is a scheme's, whose labels are names: 2 and 3 are nodes.
    path = tmp_path / "scheme.nwk"
    path.write_text("((2.1,2.2)2,(3.1,3.2)3)root;")
    assert set(read_hierarchy(str(path)).nodes) == {"root", "2", "3", "2.1", "2.2", "3.1", "3.2"}
    # The same text read as a genome tree takes them for support values.
    assert set(read_genome_tree(str(path)).support.values()) == {0.02, 0.03}


def test_two_internal_nodes_with_the_same_support_do_not_collide():
    # Read as names, the two '100' labels gave one node two parents.
    tree = Hierarchy.from_newick("((a,b)100,(c,d)100);", support_labels=True)
    assert len(tree.leaves) == 4 and len(tree) == 7


def test_an_iqtree_label_gives_the_last_support_value():
    tree = Hierarchy.from_newick("((a,b)80.5/99,c);", support_labels=True)
    assert list(tree.support.values()) == [0.99]


def test_named_internal_nodes_and_quoted_leaves_are_kept():
    tree = Hierarchy.from_newick("(('leaf one':1,'it''s':2)A.1:0.5,c)root;")
    assert {"leaf one", "it's", "A.1", "root"} <= set(tree.nodes)
    assert tree.lengths["A.1"] == 0.5


@pytest.mark.parametrize("text,message", [
    ("((a,b),);", "leaf without a name"),
    ("(a,b);(c,d);", "more than one tree"),
    ("(a,b", "ends inside a group"),
    ("(a:x,b);", "not a branch length"),
])
def test_a_malformed_newick_tree_is_refused(text, message):
    with pytest.raises(HierarchyError, match=message):
        Hierarchy.from_newick(text)


def test_repgenr2_tables_give_support_and_lengths(tmp_path):
    edges = tmp_path / "tree2tax.tsv"
    edges.write_text("child\tparent\nn1\troot\nn2\troot\nG_s_GCF_000000001.1\tn1\n"
                     "G_s_GCF_000000002.1\tn1\nG_s_GCF_000000003.1\tn2\n"
                     "G_s_GCF_000000004.1\tn2\n")
    nodes = tmp_path / "tree2tax_nodes.tsv"
    nodes.write_text("node\tparent\tn_leaves\tdepth\tbranch_length\tsupport\tred\tlineage\n"
                     "root\t\t4\t0\tNA\tNA\t0\t\nn1\troot\t2\t0.1\t0.1\t0.95\t0.4\t\n"
                     "n2\troot\t2\t0.2\t0.2\tNA\t0.5\t\n")
    tree = read_genome_tree(str(edges), str(nodes))
    assert tree.support == {"n1": 0.95}
    assert tree.lengths == {"n1": 0.1, "n2": 0.2}


# -- matching names ------------------------------------------------------------


def test_an_accession_finds_the_record_that_holds_it():
    records = ["Fam_Gen_sp1_GCF_000000001.1", "Fam.Gen.sp2.GCF_000000001.10"]
    matched, unmatched, ambiguous = resolve_names(
        ["GCF_000000001.1", "GCF_000000001.10", "Fam_Gen_sp1_GCF_000000001.1", "other"],
        records)
    assert matched["GCF_000000001.1"] == records[0], "version 1 is not version 10"
    assert matched["GCF_000000001.10"] == records[1]
    assert matched["Fam_Gen_sp1_GCF_000000001.1"] == records[0]
    assert unmatched == ["other"] and ambiguous == []


def test_an_accession_in_two_records_matches_neither():
    _matched, _unmatched, ambiguous = resolve_names(
        ["GCF_000000001.1"], ["a_GCF_000000001.1", "b_GCF_000000001.1"])
    assert ambiguous == ["GCF_000000001.1"]


# -- collapsing ----------------------------------------------------------------


def collapsed_tree():
    return Hierarchy.from_newick(
        "(((a:1,b:1)0.99:5,(c:1,d:1)0.40:5)0.95:1,((e:1,f:1)0.99:0.5,g:1)0.98:3);",
        support_labels=True)


def test_weakly_supported_nodes_are_removed_and_their_genomes_move_up():
    result = collapse(collapsed_tree(), min_support=0.9)
    assert len(result.removed) == 1
    weak = next(iter(result.removed))
    assert "support 0.4" in result.removed[weak]
    # c and d now sit at the node that held their removed clade.
    assert result.genomes["c"] == result.genomes["d"]
    assert result.genomes["c"] == result.hierarchy.parents[result.genomes["a"]]


def test_short_branches_and_small_clades_are_removed():
    by_length = collapse(collapsed_tree(), min_length=1.0)
    assert any("branch length 0.5" in reason for reason in by_length.removed.values())
    assert by_length.genomes["e"] == by_length.genomes["g"]

    by_size = collapse(collapsed_tree(), min_genomes=3)
    assert by_size.genomes["a"] == by_size.genomes["d"] != by_size.genomes["e"]


def test_folded_genomes_count_and_are_placed_with_their_leaf():
    members = {"a": ["a2", "a3"]}
    result = collapse(collapsed_tree(), min_genomes=4, members=members)
    assert result.genomes["a2"] == result.genomes["a3"] == result.genomes["a"]
    # With its two members the clade of a and b holds four genomes and stays.
    assert result.genomes["a"] != result.genomes["c"]


def test_a_node_with_the_genomes_of_its_parent_is_joined_into_it():
    tree = Hierarchy.from_newick("((((a,b)x)y,c)z);")
    result = collapse(tree)
    assert "x" not in result.hierarchy or "y" not in result.hierarchy
    assert any("same genomes as its parent" in reason for reason in result.removed.values())


def test_an_outgroup_is_placed_nowhere_and_the_root_is_the_ingroups():
    result = collapse(Hierarchy.from_newick(PLANTED_TREE, support_labels=True), outgroups=["ref"])
    assert "ref" not in result.genomes
    assert len(result.hierarchy.children[result.hierarchy.root]) == 2
    with pytest.raises(HierarchyError, match="no leaf named"):
        collapse(Hierarchy.from_newick(PLANTED_TREE, support_labels=True), outgroups=["nothere"])


def test_renaming_gives_hierarchical_names_largest_first():
    tree = Hierarchy.from_newick(PLANTED_TREE, support_labels=True)
    result = rename(collapse(tree, outgroups=["ref"]), "T")
    assert result.genomes["s1"] == "T.1.1", "A holds three genomes, B two"
    assert result.genomes["s3"] == "T.2"
    assert result.source["T"] in tree


def test_the_written_files_are_read_back_by_extract(tmp_path):
    tree = Hierarchy.from_newick(PLANTED_TREE, support_labels=True)
    result = rename(collapse(tree, outgroups=["ref"]), "T")
    hierarchy_path, genomes_path, nodes_path = write_collapsed(result, tree, str(tmp_path), "t")
    assert read_hierarchy(hierarchy_path).root == "T"
    clusters = read_genome_map(genomes_path)
    assert clusters["T.1.1"] == ["s1", "s1b"]
    assert "0.97" in Path(nodes_path).read_text()


def test_a_scheme_is_extracted_from_a_collapsed_tree(dataset, tmp_path):
    # The new scheme route: tree, collapse, extract. With the reference as the
    # outgroup, every node has a genome outside it, so even the root's canSNP
    # is found and polarised.
    tree_file = tmp_path / "tree.nwk"
    tree_file.write_text(PLANTED_TREE)
    assert cli_run("collapse", "--tree", str(tree_file), "--outgroup", "ref",
                   "--rename", "T", "--outdir", str(tmp_path), "--prefix", "t") == 0
    snps = tmp_path / "snps.vcf"
    assert cli_run("extract", "--genome-map", str(tmp_path / "t_genomes.tsv"),
                   "--hierarchy", str(tmp_path / "t_hierarchy.tsv"),
                   "--alignment", dataset.alignment, "--reference", dataset.reference,
                   "-o", str(snps)) == 0
    found = {record.node: record.pos for record in vcf_records(str(snps)) if record.node}
    # T is the planted root, T.1 is A, T.1.1 is A.1, T.2 is B.
    assert found["T"] == PLANTED["root"][0] + 1
    assert found["T.1"] == PLANTED["A"][0] + 1
    assert found["T.1.1"] == PLANTED["A.1"][0] + 1
    assert found["T.2"] == PLANTED["B"][0] + 1


def test_complementary_nodes_are_reported_as_a_missing_outgroup(caplog):
    alignment = {"a": "AAC", "b": "AAC", "c": "GAC"}
    with caplog.at_level(logging.WARNING):
        canSNP_extractor(alignment, {"L": ["a", "b"], "R": ["c"]}, snp_extractor(alignment))
    assert "whose isolates together are every genome" in caplog.text
    assert "same isolates" not in caplog.text


# -- placing genomes on an existing scheme -----------------------------------


def test_the_planted_genomes_are_placed_where_they_belong(written, dataset):
    outdir, _scores = written
    placement = place(read_scheme(outdir), Hierarchy.from_newick(PLANTED_TREE, support_labels=True),
                      read_alignment(dataset.alignment), reference_parser(dataset.reference))
    assert {node: sorted(genomes) for node, genomes in placement.clusters.items()} == {
        node: sorted(genomes) for node, genomes in CLUSTERS.items()}
    assert all(match.jaccard == 1.0 for match in placement.matches)
    assert "ref" not in {genome for genomes in placement.clusters.values() for genome in genomes}


def test_a_tree_rooted_inside_a_clade_still_places_the_genomes(written, dataset):
    # A midpoint root can fall inside a clade of the scheme. Here the tree is
    # rooted between A.1 and the rest, so A is no longer a clade of the tree
    # as rooted, only the other side of a branch.
    outdir, _scores = written
    rerooted = "((s1,s1b),(s2,((s3,s4),ref)));"
    placement = place(read_scheme(outdir), Hierarchy.from_newick(rerooted),
                      read_alignment(dataset.alignment), reference_parser(dataset.reference))
    assert sorted(placement.clusters["A.2"]) == ["s2"]
    assert sorted(placement.clusters["A.1"]) == ["s1", "s1b"]
    outside = [match for match in placement.matches if match.outside]
    assert outside and any(finding.category == "tree-root" for finding in placement.findings)


def five_and_two():
    """A clade X of five genomes and a clade Y of two, with the reference r."""
    paths = [["root"], ["root", "X"], ["root", "Y"]]
    hierarchy, _conflicts, _problems = hierarchy_from_paths(paths)
    scheme = Scheme(organism="t", paths=paths, hierarchy=hierarchy, snps=[
        SchemeSnp(node="X", strain="t", reference="t", genome="REF", position=2,
                  derived="T", ancestral="C"),
    ])
    reference = Reference.from_records([("c1", "ACGTA")])
    alignment = {"r": "ACGTA", "y1": "ACGTA", "y2": "ACGTA",
                 "x1": "ATGTA", "x2": "ATGTA", "x3": "ATGTA", "x4": "ATGTA",
                 # x5 sits in the clade but lacks the derived base
                 "x5": "ACGTA"}
    tree = Hierarchy.from_newick("(r,((x1,x2,x3,x4,x5)X,(y1,y2)Y));")
    return scheme, reference, alignment, tree


def test_a_genome_lacking_the_derived_base_stays_in_its_clade():
    # Typing would place x5 outside X, because it lacks X's derived base; the
    # tree keeps it inside, where the audit can see that it lacks it.
    scheme, reference, alignment, tree = five_and_two()
    placement = place(scheme, tree, alignment, reference)
    assert sorted(placement.clusters["X"]) == ["x1", "x2", "x3", "x4", "x5"]
    x = next(match for match in placement.matches if match.node == "X")
    assert x.jaccard == pytest.approx(0.8) and x.clade_only == 1
    assert any(finding.category == "partial-match" for finding in placement.findings)

    findings = audit_scheme(scheme, reference, alignment=alignment,
                            clusters=placement.clusters)
    failing = [finding for finding in findings if finding.category == "not-canonical"]
    assert len(failing) == 1 and "x5" in failing[0].detail


def test_a_poor_match_places_nothing():
    scheme, reference, alignment, tree = five_and_two()
    placement = place(scheme, tree, alignment, reference, min_match=0.9)
    assert "X" not in placement.clusters
    assert any(finding.category == "unplaced-node" for finding in placement.findings)


def test_genomes_folded_into_a_leaf_are_placed_with_it(written, dataset, tmp_path):
    # s1b stands in the alignment, but the tree has only s1, into which
    # dereplication folded it.
    outdir, _scores = written
    genome_map = tmp_path / "genomes_map.tsv"
    genome_map.write_text("s1\ts1\ns1b\ts1\n")
    tree = Hierarchy.from_newick("(ref,((s1,s2),(s3,s4)));")
    members = read_members(str(genome_map), tree)
    assert members == {"s1": ["s1b"]}
    placement = place(read_scheme(outdir), tree, read_alignment(dataset.alignment),
                      reference_parser(dataset.reference), members=members)
    holding = {genome: node for node, genomes in placement.clusters.items()
               for genome in genomes}
    assert holding["s1b"] == holding["s1"]


# -- the command line --------------------------------------------------------


def test_place_and_audit_by_tree_from_the_command_line(written, dataset, tmp_path):
    outdir, _scores = written
    tree_file = tmp_path / "tree.nwk"
    tree_file.write_text(PLANTED_TREE)
    genome_map = tmp_path / "placed.tsv"
    report = tmp_path / "report.tsv"
    assert cli_run("place", "--scheme-dir", outdir, "--tree", str(tree_file),
                   "--alignment", dataset.alignment, "--reference", dataset.reference,
                   "-o", str(genome_map), "--report", str(report)) == 0
    assert "s1\tA.1" in genome_map.read_text()
    assert report.read_text().splitlines()[0].startswith("node\ttree_node\tside")

    assert cli_run("audit", "--scheme-dir", outdir, "--reference", dataset.reference,
                   "--alignment", dataset.alignment, "--tree", str(tree_file)) == 0
    assert cli_run("audit", "--scheme-dir", outdir, "--reference", dataset.reference,
                   "--alignment", dataset.alignment, "--tree", str(tree_file),
                   "--genome-map", dataset.genome_map) == 1, "one or the other"


def test_collapse_reports_a_missing_file(tmp_path):
    assert cli_run("collapse", "--tree", str(tmp_path / "none.nwk"),
                   "--outdir", str(tmp_path)) == 1



def test_a_members_file_from_another_tree_is_reported(tmp_path, caplog):
    genome_map = tmp_path / "genomes_map.tsv"
    genome_map.write_text("a2\ta\nz9\tgone\n")
    tree = Hierarchy.from_newick("((a,b),c);")
    assert read_members(str(genome_map), tree) == {"a": ["a2"]}
    assert "1 leaves the tree does not have" in caplog.text and "tree2tax" in caplog.text


# -- placing genomes by the clades of the scheme's own genomes ---------------


def read_genome_tree_text(text: str) -> Hierarchy:
    """A genome tree from Newick text, as read_genome_tree reads a file."""
    return Hierarchy.from_newick(text, support_labels=True)


CLADE_PARENTS = {"A": None, "A.1": "A", "A.1.1": "A.1", "A.2": "A"}


#: Two scheme genomes at each of A.1.1 and A.2, so that both have a clade
#: wider than one leaf.
CLADE_GENOMES = {"r1": "A.1.1", "r1b": "A.1.1", "r2": "A.1", "r3": "A.2", "r3b": "A.2"}


def test_place_by_clade_deepest_node():
    tree = read_genome_tree_text("(((((r1,h1),r1b),r2),((r3,h2),r3b)),out);")
    placed, unplaced = place_by_clade(tree, CLADE_GENOMES, CLADE_PARENTS, ["out"])
    assert placed == {"h1": "A.1.1", "h2": "A.2"} and unplaced == []


def test_place_by_clade_a_node_with_one_genome_is_its_leaf_alone():
    # The clade of one genome is its leaf, as the most recent common ancestor
    # of one leaf is, and as dendropy's mrca gives it in alltree.py: h1, the
    # sister of A.1.1's only genome, may lie inside A.1.1 or beside it, and is
    # placed at A.1, the deepest node that certainly holds it.
    tree = read_genome_tree_text("((((r1,h1),r2),(r3,h2)),out);")
    placed, unplaced = place_by_clade(tree, {"r1": "A.1.1", "r2": "A.1", "r3": "A.2"},
                                      CLADE_PARENTS, ["out"])
    assert placed == {"h1": "A.1", "h2": "A"} and unplaced == []


def test_place_by_clade_roots_the_tree_on_the_outgroup():
    # Written rooted inside the ingroup; rooted on out, A.1.1 is ((r1, h1), r1b).
    tree = read_genome_tree_text("(((r1,h1),r1b),(r2,(((r3,h2),r3b),out)));")
    placed, unplaced = place_by_clade(tree, CLADE_GENOMES, CLADE_PARENTS, ["out"])
    assert placed == {"h1": "A.1.1", "h2": "A.2"} and unplaced == []


def test_place_by_clade_node_without_genomes_in_tree():
    # A.1.1's only genome r1 is missing from the tree: h1 goes to A.1
    tree = read_genome_tree_text("((((h1,r2),r2b),((r3,h2),r3b)),out);")
    placed, unplaced = place_by_clade(
        tree, {"r1": "A.1.1", "r2": "A.1", "r2b": "A.1", "r3": "A.2", "r3b": "A.2"},
        CLADE_PARENTS, ["out"])
    assert placed["h1"] == "A.1"
    assert placed["h2"] == "A.2" and unplaced == []


def test_place_by_clade_leaf_outside_every_clade():
    # h3 sits beside the outgroup
    tree = read_genome_tree_text("((((((r1,h1),r1b),r2),((r3,h2),r3b)),h3),out);")
    placed, unplaced = place_by_clade(tree, CLADE_GENOMES, CLADE_PARENTS, ["out"])
    assert "h3" in unplaced and "h3" not in placed
    assert placed == {"h1": "A.1.1", "h2": "A.2"}


def test_place_by_clade_matches_scheme_genomes_by_accession():
    tree = read_genome_tree_text(
        "((((Ec_GCF_000000001.1,h1),Ec_GCF_000000002.1),(Ec_GCF_000000003.1,"
        "Ec_GCF_000000004.1)),out);")
    placed, _unplaced = place_by_clade(
        tree, {"GCF_000000001.1": "A.1", "GCF_000000002.1": "A.1", "GCF_000000003.1": "A.2",
               "GCF_000000004.1": "A.2"}, CLADE_PARENTS, ["out"])
    assert placed == {"h1": "A.1"}


def test_place_by_clade_refuses_a_node_the_scheme_lacks():
    tree = read_genome_tree_text("((r1,h1),out);")
    with pytest.raises(ValueError, match=r"A\.9"):
        place_by_clade(tree, {"r1": "A.9"}, CLADE_PARENTS, ["out"])


def clade_scheme(tmp_path):
    """A scheme directory whose tree file is A, A.1, A.1.1 and A.2."""
    scheme_dir = tmp_path / "clade_scheme"
    scheme_dir.mkdir()
    (scheme_dir / "t_tree.txt").write_text("A\nA\tA.1\nA\tA.1\tA.1.1\nA\tA.2\n")
    (scheme_dir / "t_snp.txt").write_text(
        "snp_id\tstrain\treference\tgenome\tposition\tderived_base\tancestral_base\n"
        "A.1\tt\tt\tREF\t2\tT\tC\n")
    (scheme_dir / "t_references.txt").write_text(
        "genome\tstrain\tgenbank_id\trefseq_id\tassembly_name\nREF\tt\tGCA_1\tGCF_1\tASM1\n")
    return str(scheme_dir)


def test_place_by_clade_from_the_command_line(tmp_path):
    tree_file = tmp_path / "tree.nwk"
    tree_file.write_text("((((((r1,h1),r1b),r2),((r3,h2),r3b)),h3),out);")
    scheme_map = tmp_path / "scheme_genomes.tsv"
    scheme_map.write_text("".join(f"{genome}\t{node}\n" for genome, node in CLADE_GENOMES.items()))
    # h4 is a genome of the alignment that the tree does not hold.
    alignment = tmp_path / "sites.fasta"
    alignment.write_text("".join(f">{name}\nA\n" for name in
                                 (*CLADE_GENOMES, "h1", "h2", "h3", "h4", "out")))
    output = tmp_path / "placed.tsv"
    report = tmp_path / "report.tsv"
    assert cli_run("place", "--rule", "clade", "--scheme-dir", clade_scheme(tmp_path),
                   "--tree", str(tree_file), "--genome-map", str(scheme_map),
                   "--outgroup", "out", "--alignment", str(alignment),
                   "-o", str(output), "--report", str(report)) == 0
    placed = read_genome_map(str(output))
    assert {node: sorted(genomes) for node, genomes in placed.items()} == {
        "A.1.1": ["h1", "r1", "r1b"], "A.1": ["r2"], "A.2": ["h2", "r3", "r3b"]}
    rows = [line.split("\t") for line in report.read_text().splitlines()]
    assert rows[0] == ["genome", "node", "note"]
    notes = {row[0]: (row[1], row[2]) for row in rows[1:]}
    assert notes["h1"] == ("A.1.1", "placed by the tree")
    assert notes["r1"] == ("A.1.1", "scheme genome")
    assert notes["h3"][0] == "" and "no clade" in notes["h3"][1]
    assert notes["h4"][0] == "" and "not a leaf" in notes["h4"][1]
    assert "out" not in notes


def test_place_by_clade_needs_a_genome_map_and_an_outgroup(tmp_path):
    tree_file = tmp_path / "tree.nwk"
    tree_file.write_text("((r1,h1),out);")
    scheme_dir = clade_scheme(tmp_path)
    with pytest.raises(SystemExit):
        cli_run("place", "--rule", "clade", "--scheme-dir", scheme_dir, "--tree", str(tree_file),
                "--outgroup", "out", "-o", str(tmp_path / "out.tsv"))
    with pytest.raises(SystemExit):
        cli_run("place", "--rule", "clade", "--scheme-dir", scheme_dir, "--tree", str(tree_file),
                "--genome-map", str(tree_file), "-o", str(tmp_path / "out.tsv"))
    with pytest.raises(SystemExit):
        cli_run("place", "--scheme-dir", scheme_dir, "--tree", str(tree_file),
                "-o", str(tmp_path / "out.tsv"))


def test_audit_by_clade_from_the_command_line(written, dataset, tmp_path):
    outdir, _scores = written
    tree_file = tmp_path / "tree.nwk"
    tree_file.write_text(PLANTED_TREE)
    scheme_map = tmp_path / "scheme_genomes.tsv"
    scheme_map.write_text("s1\tA.1\ns1b\tA.1\ns2\tA.2\ns3\tB.1\ns4\tB.2\n")
    clade = ["--tree", str(tree_file), "--rule", "clade", "--outgroup", "ref"]
    common = ["audit", "--scheme-dir", outdir, "--reference", dataset.reference,
              "--alignment", dataset.alignment]
    assert cli_run(*common, *clade, "--genome-map", str(scheme_map)) == 0

    # Without s1b among the scheme's genomes, A.1's clade is s1's leaf alone,
    # and the tree places s1b at A, where the audit finds A.1's derived base.
    scheme_map.write_text("s1\tA.1\ns2\tA.2\ns3\tB.1\ns4\tB.2\n")
    findings = tmp_path / "findings.tsv"
    assert cli_run(*common, *clade, "--genome-map", str(scheme_map), "-o", str(findings)) == 2
    assert "s1b; they are placed at A " in findings.read_text()

    assert cli_run(*common, "--genome-map", str(scheme_map), "--rule", "clade",
                   "--outgroup", "ref") == 1, "clade needs --tree"
    assert cli_run(*common, "--tree", str(tree_file), "--outgroup", "ref") == 1, \
        "--outgroup applies to clade"


def test_place_by_clade_ties_in_depth_go_to_the_smaller_clade():
    # B1 and B2 are siblings. B1's genomes enclose B2's, so h lies in both
    # clades; B2's is the smaller (3 leaves against 5), and B1 would come
    # first by name.
    tree = read_genome_tree_text("(((b1a,((b2a,h),b2b)),b1b),out);")
    parents = {"A": None, "B1": "A", "B2": "A"}
    placed, unplaced = place_by_clade(
        tree, {"b1a": "B1", "b1b": "B1", "b2a": "B2", "b2b": "B2"}, parents, ["out"])
    assert placed == {"h": "B2"} and unplaced == []


def test_place_by_clade_warns_of_a_scheme_genome_that_is_not_a_leaf(tmp_path, caplog):
    from cansnpscorer.cansnp_tree import place_clade

    tree_file = tmp_path / "tree.nwk"
    tree_file.write_text("(((((r1,h1),r1b),r2),((r3,h2),r3b)),out);")
    scheme_map = tmp_path / "scheme_genomes.tsv"
    genomes = {**CLADE_GENOMES, "r4": "A.2"}
    scheme_map.write_text("".join(f"{genome}\t{node}\n" for genome, node in genomes.items()))
    scheme = read_scheme(clade_scheme(tmp_path))
    with caplog.at_level(logging.WARNING):
        placement = place_clade(scheme, str(tree_file), str(scheme_map), ["out"],
                                [*genomes, "h1", "h2", "out"])
    assert "1 scheme genomes of the alignment are not leaves of the tree" in caplog.text
    assert "r4" in caplog.text
    assert "r4" not in {g for genomes in placement.clusters.values() for g in genomes}


def test_audit_by_clade_does_not_flag_a_reference_the_tree_lacks(written, dataset, tmp_path,
                                                                caplog):
    outdir, _scores = written
    tree_file = tmp_path / "tree.nwk"
    tree_file.write_text("(((s1:0.01,s1b:0.01):0.2,s2:0.3):0.4,(s3:0.2,s4:0.2):0.4,out:0.5);")
    scheme_map = tmp_path / "scheme_genomes.tsv"
    scheme_map.write_text("s1\tA.1\ns1b\tA.1\ns2\tA.2\ns3\tB.1\ns4\tB.2\n")
    with caplog.at_level(logging.WARNING):
        cli_run("audit", "--scheme-dir", outdir, "--reference", dataset.reference,
                "--alignment", dataset.alignment, "--tree", str(tree_file), "--rule", "clade",
                "--outgroup", "out", "--genome-map", str(scheme_map))
    assert "not leaves of the tree" not in caplog.text


def test_no_new_clades_needs_an_alignment(written, dataset):
    outdir, _scores = written
    assert cli_run("audit", "--scheme-dir", outdir, "--reference", dataset.reference,
                   "--no-new-clades") == 1


def test_place_by_clade_refuses_a_min_match(tmp_path):
    tree_file = tmp_path / "tree.nwk"
    tree_file.write_text("((r1,h1),out);")
    with pytest.raises(SystemExit):
        cli_run("place", "--rule", "clade", "--scheme-dir", clade_scheme(tmp_path),
                "--tree", str(tree_file), "--genome-map", str(tree_file), "--outgroup", "out",
                "--min-match", "0.8", "-o", str(tmp_path / "out.tsv"))
