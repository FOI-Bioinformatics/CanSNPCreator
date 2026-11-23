"""Defects of collapse and place found in the fourth audit round, each pinned
by the input that showed it."""

import time

import pytest

from cansnpscorer.hierarchy import Hierarchy, HierarchyError
from cansnpscorer.model import Reference
from cansnpscorer.scheme import Scheme, SchemeSnp, hierarchy_from_paths
from cansnpscorer.trees.collapsing import collapse, rename
from cansnpscorer.trees.placing import place
from cansnpscorer.trees.reading import read_genome_tree, read_members
from cansnpscorer.trees.rooting import reroot


def newick(text: str) -> Hierarchy:
    return Hierarchy.from_newick(text, support_labels=True)


def caterpillar(depth: int) -> str:
    """A tree in which every internal node has one leaf and one internal child."""
    return "(" * depth + "L0" + "".join(f",L{i}):0.1" for i in range(1, depth + 1)) + ";"


# -- rooting -----------------------------------------------------------------


def test_a_root_with_one_child_is_not_a_genome_after_rerooting():
    # The outer brackets leave the root with one child. Rooted on O, the old
    # root became a leaf and was placed as a genome.
    tree = newick("(((A:1,O:1)70:1,(B:1,C:1)90:2):7);")
    rooted = reroot(tree, ["O"])
    assert sorted(rooted.leaves) == ["A", "B", "C", "O"]
    collapsed = collapse(tree, outgroups=["O"])
    assert sorted(collapsed.genomes) == ["A", "B", "C"]
    assert collapsed.genomes["B"] == collapsed.genomes["C"] != collapsed.genomes["A"]


# -- reading Newick ------------------------------------------------------------


def test_a_deep_tree_is_read_collapsed_and_used_for_placing_without_recursion():
    # A caterpillar of 3000 genomes exceeded Python's recursion limit.
    tree = newick(caterpillar(3000))
    assert len(tree.leaves) == 3001
    collapsed = rename(collapse(tree, outgroups=["L0"]), "D")
    assert len(collapsed.genomes) == 3000 and "L0" not in collapsed.genomes

    paths = [["root"], ["root", "X"]]
    hierarchy, _, _ = hierarchy_from_paths(paths)
    scheme = Scheme(
        organism="t",
        paths=paths,
        hierarchy=hierarchy,
        snps=[
            SchemeSnp(
                node="X",
                strain="t",
                reference="t",
                genome="REF",
                sequence="c",
                position=2,
                derived="G",
                ancestral="A",
            )
        ],
    )
    # The 10 genomes deepest in the caterpillar carry the derived base.
    alignment = {"ref": "AAA"}
    alignment.update({f"L{i}": "AGA" if i <= 10 else "AAA" for i in range(3001)})
    placement = place(scheme, tree, alignment, Reference.from_records([("c", "AAA")]))
    assert sorted(placement.clusters["X"]) == sorted(f"L{i}" for i in range(11))


def test_a_newick_tree_opening_with_a_comment_is_read(tmp_path):
    # Refused as neither Newick nor a table, since it does not start with '('.
    path = tmp_path / "tree.nwk"
    path.write_text("[&R] ((A,B)90,C);\n")
    assert read_genome_tree(str(path)).support == {"node1": 0.9}
    path.write_text("[a comment\nover two lines]\n((A,B)90,C);\n")
    assert sorted(read_genome_tree(str(path)).leaves) == ["A", "B", "C"]


@pytest.mark.parametrize(
    "label, value",
    [
        ("1e-1", 0.1),
        ("1.", 1.0),
        (".95", 0.95),
        ("/90", 0.9),
        ("85.3/97", 0.97),
        ("0.95", 0.95),
        ("95/100", 1.0),
    ],
)
def test_support_labels_in_every_number_form_are_support(label, value):
    # Read as names, these made a support value a node of the scheme.
    tree = newick(f"((A,B){label},C);")
    assert label not in tree
    assert tree.support == {"node1": pytest.approx(value)}


def test_a_support_label_without_its_last_value_is_no_name_and_no_support():
    tree = newick("((A,B)80/,C);")
    assert "80/" not in tree and tree.support == {}


# -- renaming ------------------------------------------------------------------


def test_rename_orders_children_of_equal_size_by_their_genomes():
    # The names followed the order of the branches in the tree file.
    first = "(O:1,((A:1,B:1)95:1,(C:1,D:1)96:1,(E:1,F:1,G:1)97:1)98:1);"
    second = "(O:1,((C:1,D:1)96:1,(A:1,B:1)95:1,(E:1,F:1,G:1)97:1)98:1);"
    named = [rename(collapse(newick(text), outgroups=["O"]), "T").genomes for text in (first, second)]
    assert named[0] == named[1]
    assert named[0]["E"] == "T.1", "largest first"
    assert named[0]["A"] == "T.2" and named[0]["C"] == "T.3", "then by the first genome"


# -- outgroups -----------------------------------------------------------------


@pytest.mark.parametrize("outgroup", ["GCF_000007125", "GCF_000007125.1", "Brucella_sp_GCF_000007125.1"])
def test_an_outgroup_given_by_accession_matches_its_leaf(outgroup):
    # Only the full leaf name was accepted.
    tree = newick("(Brucella_sp_GCF_000007125.1:1,((A:1,B:1)95:1,(C:1,D:1)96:1)98:1);")
    collapsed = collapse(tree, outgroups=[outgroup])
    assert sorted(collapsed.genomes) == ["A", "B", "C", "D"]
    assert reroot(tree, [outgroup]) is not None


def test_an_outgroup_matching_two_leaves_is_refused():
    tree = newick("(x_GCF_000007125.1,(y_GCF_000007125.2,(A,B)95)98);")
    with pytest.raises(HierarchyError, match="more than one leaf"):
        collapse(tree, outgroups=["GCF_000007125"])


# -- scaling -------------------------------------------------------------------


def test_collapse_place_and_read_members_take_linear_time(tmp_path):
    # Each was quadratic: 8000 leaves took 19 s to collapse, and 20000 rows
    # of members 100 s to read. The bounds are generous, so that only a
    # return to quadratic time fails them.
    start = time.perf_counter()
    tree = newick(caterpillar(20000))
    collapse(tree, min_genomes=2, outgroups=["L0"])
    leaves = [f"Fam_Gen_sp_GCF_{i:09d}.1" for i in range(2000)]
    flat = Hierarchy([(leaf, "r") for leaf in leaves])
    path = tmp_path / "genomes_map.tsv"
    path.write_text("".join(f"GCF_{100000 + i:09d}.1\t{leaves[i % 2000]}\n" for i in range(20000)))
    members = read_members(str(path), flat)
    assert sum(len(names) for names in members.values()) == 20000
    assert time.perf_counter() - start < 20


def test_a_run_of_single_children_is_joined_into_its_highest_node():
    # Joining is now done in one pass; the result is the same as joining one
    # node at a time.
    collapsed = collapse(newick("(O:1,((((A:1,B:1)95:1)96:1)97:1,C:1)99:1);"), outgroups=["O"])
    assert collapsed.genomes == {"A": "node3", "B": "node3", "C": "node5"}
    # The root holds no genome and node4 is its one child, the outgroup
    # being at no node.
    assert list(collapsed.removed) == ["node4", "node2", "node1"]
