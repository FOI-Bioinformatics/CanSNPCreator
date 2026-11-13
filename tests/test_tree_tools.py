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
    records = ["FamGensp1GCF0000000011", "FamGensp2GCF00000000110"]
    matched, unmatched, ambiguous = resolve_names(
        ["GCF_000000001.1", "GCF_000000001.10", "Fam_Gen_sp1_GCF_000000001.1", "other"],
        records)
    assert matched["GCF_000000001.1"] == records[0], "version 1 is not version 10"
    assert matched["GCF_000000001.10"] == records[1]
    assert matched["Fam_Gen_sp1_GCF_000000001.1"] == records[0]
    assert unmatched == ["other"] and ambiguous == []


def test_an_accession_in_two_records_matches_neither():
    _matched, _unmatched, ambiguous = resolve_names(
        ["GCF_000000001.1"], ["aGCF0000000011", "bGCF0000000011"])
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
