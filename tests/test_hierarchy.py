"""The scheme hierarchy, in each of the shapes it arrives in.

The defect these guard against: the child and parent table was handed to the
walker written for the nested form, which read only its keys. Every internal
node came out with no isolates, so no canSNP could be assigned to a branch,
which is the whole purpose of the tool.
"""

import pytest
from conftest import CLUSTERS, PLANTED, node_isolates

from cansnpscorer.helpers.data_utils import canSNP_extractor, snp_extractor
from cansnpscorer.helpers.file_utils import (
    HCParser,
    read_alignment,
    read_hc_file,
    recognize_group_file_format,
    recognize_hc_file_format,
)
from cansnpscorer.hierarchy import Hierarchy, HierarchyError

NESTED = [{"root": [{"A": ["A.1", "A.2"]}, {"B": ["B.1", "B.2"]}]}]
TABLE = [("A", "root"), ("B", "root"), ("A.1", "A"), ("A.2", "A"), ("B.1", "B"), ("B.2", "B")]
NEWICK = "((A.1,A.2)A,(B.1,B.2)B)root;"

NODES = ["root", "A", "A.1", "A.2", "B", "B.1", "B.2"]


def shape(hierarchy):
    return {node: sorted(hierarchy.children[node]) for node in hierarchy.nodes}


# -- the three shapes agree --------------------------------------------------


@pytest.mark.parametrize("built", [
    pytest.param(Hierarchy.from_nested(NESTED), id="nested"),
    pytest.param(Hierarchy.from_child_parent(TABLE), id="child_parent"),
    pytest.param(Hierarchy.from_newick(NEWICK), id="newick"),
])
def test_every_shape_gives_the_same_tree(built):
    assert built.root == "root"
    assert sorted(built.nodes) == sorted(NODES)
    assert shape(built) == {
        "root": ["A", "B"], "A": ["A.1", "A.2"], "B": ["B.1", "B.2"],
        "A.1": [], "A.2": [], "B.1": [], "B.2": [],
    }
    assert built.leaves == ["A.1", "A.2", "B.1", "B.2"]
    assert built.path("B.2") == ["root", "B", "B.2"]


@pytest.mark.parametrize("built", [
    pytest.param(Hierarchy.from_nested(NESTED), id="nested"),
    pytest.param(Hierarchy.from_child_parent(TABLE), id="child_parent"),
])
def test_an_internal_node_holds_the_isolates_below_it(built):
    isolates = built.isolates(CLUSTERS)
    for node in NODES:
        assert sorted(isolates[node]) == sorted(node_isolates(node)), node
    assert sorted(isolates["root"]) == ["s1", "s1b", "s2", "s3", "s4"]


def test_a_leaf_without_isolates_is_reported_as_empty():
    built = Hierarchy.from_child_parent(TABLE)
    isolates = built.isolates({"A.1": ["s1"]})
    assert isolates["A.2"] == []
    assert isolates["A"] == ["s1"]


def test_a_cluster_that_is_not_in_the_hierarchy_is_refused():
    built = Hierarchy.from_child_parent(TABLE)
    with pytest.raises(HierarchyError, match="not in the hierarchy: C.9"):
        built.isolates({"A.1": ["s1"], "C.9": ["s9"]})


# -- what is not a tree ------------------------------------------------------


def test_a_node_with_two_parents_is_refused():
    with pytest.raises(HierarchyError, match="more than one parent"):
        Hierarchy.from_child_parent([("A", "root"), ("B", "root"), ("X", "A"), ("X", "B")])


def test_a_cycle_is_refused():
    with pytest.raises(HierarchyError, match="no root|not a tree"):
        Hierarchy.from_child_parent([("A", "B"), ("B", "A")])


def test_two_roots_are_refused():
    with pytest.raises(HierarchyError, match="2 roots where one is expected"):
        Hierarchy.from_nested([{"one": ["a"]}, {"two": ["b"]}])


def test_a_node_that_is_its_own_parent_is_refused():
    with pytest.raises(HierarchyError, match="its own parent"):
        Hierarchy.from_child_parent([("A", "A")])


def test_an_empty_hierarchy_is_refused():
    with pytest.raises(HierarchyError, match="holds no nodes"):
        Hierarchy.from_nested([])


# -- reading the files -------------------------------------------------------


def test_both_hierarchy_files_of_the_dataset_describe_the_same_tree(dataset):
    assert recognize_hc_file_format(dataset.hc_yaml) == "yaml"
    assert recognize_hc_file_format(dataset.hc_repgenr) == "repgenr"
    from_yaml = read_hc_file(dataset.hc_yaml)
    from_table = read_hc_file(dataset.hc_repgenr)
    assert shape(from_yaml) == shape(from_table)


def test_the_repgenr_route_finds_the_same_cansnps_as_the_nested_one(dataset):
    # The defect itself: through the repgenr files, every canSNP on an internal
    # branch used to be missed.
    alignment = read_alignment(dataset.alignment)
    snps = snp_extractor(alignment)

    found = {}
    for hc_file, group_file in [(dataset.hc_yaml, dataset.groups_csv),
                                (dataset.hc_repgenr, dataset.groups_repgenr)]:
        hierarchy = read_hc_file(hc_file)
        clusters = HCParser(group_file, data_format=recognize_group_file_format(group_file))
        cansnps = canSNP_extractor(alignment, hierarchy.isolates(clusters.isolates), snps)
        found[hc_file] = {cansnp.column: cansnp.node for cansnp in cansnps if cansnp.node}

    through_yaml, through_table = found.values()
    assert through_yaml == through_table
    assert through_table == {column: node for node, (column, _a, _d) in PLANTED.items()
                             for node in [node]}
    assert "A" in through_table.values() and "root" in through_table.values(), \
        "internal branches must be marked"


def test_an_unreadable_hierarchy_file_names_the_formats(tmp_path):
    path = tmp_path / "hc.txt"
    path.write_text("just a line of text\n")
    with pytest.raises(HierarchyError, match="supported are YAML, JSON, Newick"):
        recognize_hc_file_format(str(path))


def test_an_empty_hierarchy_file_is_refused(tmp_path):
    path = tmp_path / "hc.txt"
    path.write_text("# only a comment\n\n")
    with pytest.raises(HierarchyError, match="is empty, or holds only comments"):
        recognize_hc_file_format(str(path))


def test_a_repgenr_table_without_its_columns_is_refused(tmp_path):
    path = tmp_path / "hc.tsv"
    path.write_text("child\tmother\nA\troot\n")
    with pytest.raises(HierarchyError, match="needs 'child' and 'parent' columns"):
        read_hc_file(str(path), "repgenr")


def test_a_json_hierarchy_is_read(tmp_path):
    path = tmp_path / "hc.json"
    path.write_text('{"root": {"A": ["A.1", "A.2"], "B": ["B.1", "B.2"]}}')
    assert shape(read_hc_file(str(path))) == shape(Hierarchy.from_nested(NESTED))


def test_a_newick_hierarchy_is_read(tmp_path):
    path = tmp_path / "hc.newick"
    path.write_text("((A.1:0.1,A.2:0.2)A:0.3,(B.1,B.2)B)root;\n")
    assert recognize_hc_file_format(str(path)) == "newick"
    assert shape(read_hc_file(str(path))) == shape(Hierarchy.from_nested(NESTED))


def test_an_unnamed_internal_node_is_given_a_name(tmp_path):
    path = tmp_path / "hc.newick"
    path.write_text("((A.1,A.2),(B.1,B.2))root;\n")
    built = read_hc_file(str(path))
    assert built.root == "root"
    assert sorted(built.leaves) == ["A.1", "A.2", "B.1", "B.2"]
    assert len(built) == 7


# -- group files -------------------------------------------------------------


def test_the_group_columns_may_be_in_either_order(tmp_path):
    path = tmp_path / "groups.csv"
    path.write_text("Group,Strain\nA.1,s1\nA.1,s1b\n")
    assert recognize_group_file_format(str(path)) == "csv"
    assert HCParser(str(path), data_format="csv").isolates == {"A.1": ["s1", "s1b"]}


def test_a_group_file_without_its_columns_is_refused(tmp_path):
    path = tmp_path / "groups.csv"
    path.write_text("Isolate,Cluster\ns1,A.1\n")
    with pytest.raises(ValueError, match="cannot tell which group file format"):
        recognize_group_file_format(str(path))

    with pytest.raises(ValueError, match="needs 'Strain' and 'Group' columns"):
        HCParser(str(path), data_format="csv")


def test_cluster_names_keep_their_punctuation(dataset):
    # Cluster names must read exactly as the hierarchy names them, so A.1 may
    # not be cleaned to A1, while genome names are cleaned to match the
    # alignment.
    clusters = HCParser(dataset.groups_repgenr, data_format="repgenr_group").isolates
    assert set(clusters) == set(CLUSTERS)
