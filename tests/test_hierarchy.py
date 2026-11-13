"""The scheme hierarchy, and the genome map that places genomes on it.

The hierarchy is read from a child and parent table, as FlexTaxD and RepGenR2's
tree2tax.tsv hold one, or from Newick. The genome map is two tab separated
columns, genome and node, as FlexTaxD's --genomeid2taxid. The defect these
guard against first: an internal node must hold the isolates of every node
below it, or no canSNP can be assigned to a branch.
"""

import pytest
from conftest import CLUSTERS, node_isolates

from cansnpscorer.helpers.file_utils import hierarchy_format, read_genome_map, read_hierarchy
from cansnpscorer.hierarchy import Hierarchy, HierarchyError

TABLE = [("A", "root"), ("B", "root"), ("A.1", "A"), ("A.2", "A"), ("B.1", "B"), ("B.2", "B")]
NEWICK = "((A.1,A.2)A,(B.1,B.2)B)root;"

NODES = ["root", "A", "A.1", "A.2", "B", "B.1", "B.2"]


def shape(hierarchy):
    return {node: sorted(hierarchy.children[node]) for node in hierarchy.nodes}


# -- the two formats agree ---------------------------------------------------


@pytest.mark.parametrize("built", [
    pytest.param(Hierarchy.from_child_parent(TABLE), id="child_parent"),
    pytest.param(Hierarchy.from_newick(NEWICK), id="newick"),
])
def test_both_formats_give_the_same_tree(built):
    assert built.root == "root"
    assert sorted(built.nodes) == sorted(NODES)
    assert shape(built) == {
        "root": ["A", "B"], "A": ["A.1", "A.2"], "B": ["B.1", "B.2"],
        "A.1": [], "A.2": [], "B.1": [], "B.2": [],
    }
    assert built.leaves == ["A.1", "A.2", "B.1", "B.2"]
    assert built.path("B.2") == ["root", "B", "B.2"]


def test_the_dataset_files_describe_the_same_tree(dataset):
    assert hierarchy_format(dataset.hierarchy) == "table"
    assert hierarchy_format(dataset.hierarchy_newick) == "newick"
    assert shape(read_hierarchy(dataset.hierarchy)) == shape(
        read_hierarchy(dataset.hierarchy_newick))


def test_an_internal_node_holds_the_isolates_below_it():
    isolates = Hierarchy.from_child_parent(TABLE).isolates(CLUSTERS)
    for node in NODES:
        assert sorted(isolates[node]) == sorted(node_isolates(node)), node
    assert sorted(isolates["root"]) == ["s1", "s1b", "s2", "s3", "s4"]


def test_a_leaf_without_isolates_is_reported_as_empty():
    isolates = Hierarchy.from_child_parent(TABLE).isolates({"A.1": ["s1"]})
    assert isolates["A.2"] == []
    assert isolates["A"] == ["s1"]


def test_a_cluster_that_is_not_in_the_hierarchy_is_refused():
    with pytest.raises(HierarchyError, match="not in the hierarchy: C.9"):
        Hierarchy.from_child_parent(TABLE).isolates({"A.1": ["s1"], "C.9": ["s9"]})


# -- what is not a tree ------------------------------------------------------


def test_a_node_with_two_parents_is_refused():
    with pytest.raises(HierarchyError, match="more than one parent"):
        Hierarchy.from_child_parent([("A", "root"), ("B", "root"), ("X", "A"), ("X", "B")])


def test_a_cycle_is_refused():
    with pytest.raises(HierarchyError, match="no root|not a tree"):
        Hierarchy.from_child_parent([("A", "B"), ("B", "A")])


def test_two_roots_are_refused():
    with pytest.raises(HierarchyError, match="2 roots where one is expected"):
        Hierarchy.from_child_parent([("a", "one"), ("b", "two")])


def test_a_node_that_is_its_own_parent_is_refused():
    with pytest.raises(HierarchyError, match="its own parent"):
        Hierarchy.from_child_parent([("A", "A")])


# -- the table ---------------------------------------------------------------


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


def test_numeric_node_names_stay_names(tmp_path):
    # Read as YAML, 2.10 became the number 2.1 and the node vanished into 2.1.
    path = write(tmp_path, "h.tsv", "child\tparent\n2\troot\n2.1\t2\n2.10\t2\n1e3\troot\n")
    assert set(read_hierarchy(path).nodes) == {"root", "2", "2.1", "2.10", "1e3"}


def test_a_root_row_with_an_empty_parent_is_read(tmp_path):
    path = write(tmp_path, "h.tsv", "child\tparent\nroot\t\nA\troot\n")
    hierarchy = read_hierarchy(path)
    assert hierarchy.root == "root" and hierarchy.children["root"] == ["A"]


def test_tree2tax_columns_in_any_order_and_extra_columns_are_read(tmp_path):
    # RepGenR2 adds a rank column with --ranks.
    path = write(tmp_path, "tree2tax.tsv",
                 "parent\tchild\trank\nroot\tA\tno rank\nA\tg1\tstrain\n")
    assert read_hierarchy(path).path("g1") == ["root", "A", "g1"]


def test_windows_line_endings_are_read(tmp_path):
    path = write(tmp_path, "h.tsv", "child\tparent\r\nA\troot\r\n")
    assert read_hierarchy(path).nodes == ["root", "A"]


@pytest.mark.parametrize("text,message", [
    ("", "is empty"),
    ("- root:\n  - A\n", "neither a Newick tree nor"),
    ('{"root": ["A"]}\n', "neither a Newick tree nor"),
    ("child\tparent\n", "holds no rows"),
    ("child\tparent\n\troot\n", "names no child"),
])
def test_anything_else_is_refused(tmp_path, text, message):
    with pytest.raises(HierarchyError, match=message):
        read_hierarchy(write(tmp_path, "h.txt", text))


# -- Newick --------------------------------------------------------------------


def test_an_unnamed_internal_node_is_given_a_name(tmp_path):
    built = read_hierarchy(write(tmp_path, "h.nwk", "((A.1,A.2),(B.1,B.2)B)root;"))
    assert "node1" in built
    assert built.children["node1"] == ["A.1", "A.2"]


# -- the genome map ------------------------------------------------------------


def test_the_genome_map_places_genomes_at_nodes(dataset):
    assert read_genome_map(dataset.genome_map) == CLUSTERS


def test_genome_names_keep_their_punctuation_and_lose_a_fasta_suffix(tmp_path):
    path = write(tmp_path, "map.tsv", "# comment\nGCF_000001.1.fna.gz\tA.1\nstrain-7\tA\n")
    assert read_genome_map(path) == {"A.1": ["GCF_000001.1"], "A": ["strain-7"]}


@pytest.mark.parametrize("text,message", [
    ("Strain,Group\ns1,A.1\n", "a genome and a node"),
    ("s1\tA.1\textra\n", "a genome and a node"),
    ("s1\tA.1\ns1\tB\n", "already placed at A.1"),
    ("# only a comment\n", "places no genome"),
])
def test_a_malformed_genome_map_is_refused(tmp_path, text, message):
    with pytest.raises(ValueError, match=message):
        read_genome_map(write(tmp_path, "map.tsv", text))
