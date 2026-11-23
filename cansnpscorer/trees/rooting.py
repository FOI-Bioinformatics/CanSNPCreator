"""Rooting a genome tree on its outgroup."""

import logging
from collections.abc import Sequence

from cansnpscorer.helpers.names import resolve_names
from cansnpscorer.hierarchy import Hierarchy, HierarchyError


def outgroup_leaves(tree: Hierarchy, outgroups: Sequence[str]) -> list[str]:
    """The leaves of the tree the outgroups name.

    An outgroup is matched to a leaf as genome names are everywhere else, by
    ``resolve_names``: an assembly accession, with or without its version,
    matches the one leaf whose name holds it, as RepGenR2's
    outgroup_accession.txt gives it.

    Raises:
        HierarchyError: If an outgroup matches no leaf, or more than one.
    """
    if not outgroups:
        return []
    matched, unmatched, ambiguous = resolve_names(list(dict.fromkeys(outgroups)), tree.leaves)
    if unmatched:
        raise HierarchyError(f"the tree has no leaf named {', '.join(sorted(unmatched))}")
    if ambiguous:
        raise HierarchyError(
            f"the outgroups {', '.join(sorted(ambiguous))} each match more than one leaf of the "
            f"tree; give the full leaf name"
        )
    return list(dict.fromkeys(matched.values()))


def _without_unary_root(tree: Hierarchy) -> Hierarchy:
    """The tree without a root that has a single internal child.

    Such a root adds a branch above the whole tree that splits no genomes;
    rooted elsewhere, it would become a leaf that is no genome.
    """
    root = tree.root
    removed = []
    while len(tree.children[root]) == 1 and tree.children[tree.children[root][0]]:
        removed.append(root)
        root = tree.children[root][0]
    if not removed:
        return tree
    gone = set(removed)
    stripped = Hierarchy(
        [(node, None if node == root else tree.parents[node]) for node in tree.nodes if node not in gone]
    )
    stripped.support = {node: value for node, value in tree.support.items() if node not in gone}
    # The branch above the new root led to the removed one.
    stripped.lengths = {
        node: value for node, value in tree.lengths.items() if node not in gone and node != root
    }
    return stripped


def reroot(tree: Hierarchy, outgroups: Sequence[str]) -> Hierarchy:
    """Root the tree on the branch that separates the outgroups from the rest.

    A tree builder such as IQ-TREE writes an unrooted tree rooted at its first
    taxon, so the outgroups need not sit at the root. The tree is returned
    unchanged when they already form one child of the root. Otherwise a new
    root is placed on the branch above the outgroups' clade, halving its
    length. The support and length of each branch stay with that branch,
    which, below the new root, is now held by the node at its lower end.

    A root with a single child, as some tree writers leave, is removed: on
    the rooted tree it would be a leaf.

    Raises:
        HierarchyError: If an outgroup matches no leaf of the tree, or the
            outgroups do not form a clade on any rooting of the tree.
    """
    outgroup = set(outgroup_leaves(tree, outgroups))
    if not outgroup:
        return tree
    tree = _without_unary_root(tree)
    nodes = tree.nodes
    root = tree.root
    # The leaves, and the outgroups, at or below each node.
    leaves_below: dict[str, int] = {}
    outgroups_below: dict[str, int] = {}
    for node in reversed(nodes):
        children = tree.children[node]
        leaves_below[node] = sum(leaves_below[child] for child in children) if children else 1
        outgroups_below[node] = sum(outgroups_below[child] for child in children) + (node in outgroup)

    n_ingroup = leaves_below[root] - len(outgroup)
    if not n_ingroup:
        raise HierarchyError("every leaf of the tree is an outgroup")

    def splits(node: str) -> bool:
        """Whether the clade below the node is the outgroups or the ingroup."""
        return (leaves_below[node] == outgroups_below[node] == len(outgroup)) or (
            outgroups_below[node] == 0 and leaves_below[node] == n_ingroup
        )

    # The branch that splits the tree into outgroups and ingroup is held by
    # the node below it, on whichever side of the split that node is.
    branch = next((node for node in nodes if node != root and splits(node)), None)
    if branch is None:
        raise HierarchyError(
            f"the outgroups ({', '.join(sorted(outgroup))}) do not form a clade of the tree, "
            f"however it is rooted, so no root separates them from the other genomes"
        )
    if tree.parents[branch] == root and (
        outgroups_below[branch] == len(outgroup) or len(tree.children[root]) == 2
    ):
        return tree

    # Every branch as an undirected edge, labelled by the node below it.
    neighbours: dict[str, list[tuple[str, str]]] = {node: [] for node in nodes}
    for node in nodes:
        parent = tree.parents[node]
        if parent is not None:
            neighbours[node].append((parent, node))
            neighbours[parent].append((node, node))

    new_root = "root"
    taken = set(nodes)
    number = 1
    while new_root in taken:
        number += 1
        new_root = f"root{number}"

    upper = tree.parents[branch]
    # The branch is not the root, so it has a parent.
    assert upper is not None
    pairs: list[tuple[str, str | None]] = [(new_root, None), (branch, new_root), (upper, new_root)]
    # The branch that held the root is split in two; both halves are the same
    # split of the genomes, so they carry the same support.
    support: dict[str, float] = {}
    lengths: dict[str, float] = {}
    if branch in tree.support:
        support[branch] = support[upper] = tree.support[branch]
    if branch in tree.lengths:
        lengths[branch] = lengths[upper] = tree.lengths[branch] / 2

    # Walk away from the new root, orienting every other branch away from it.
    stack = [(branch, upper), (upper, branch)]
    while stack:
        node, came_from = stack.pop()
        for other, edge in neighbours[node]:
            if other == came_from:
                continue
            pairs.append((other, node))
            if edge in tree.support:
                support[other] = tree.support[edge]
            if edge in tree.lengths:
                lengths[other] = tree.lengths[edge]
            stack.append((other, node))

    # The old root, when it had two children, is now a node with one child:
    # its two branches are one, and it is removed.
    new_children: dict[str, list[str]] = {}
    for child, parent in pairs:
        if parent is not None:
            new_children.setdefault(parent, []).append(child)
    old = root
    if old != upper and len(new_children.get(old, [])) == 1:
        (only,) = new_children[old]
        parent = next(parent for child, parent in pairs if child == old)
        pairs = [(child, parent if child == only else above) for child, above in pairs if child != old]
        if old in lengths or only in lengths:
            lengths[only] = lengths.pop(old, 0.0) + lengths.get(only, 0.0)
        value = support.pop(old, None)
        if value is not None:
            support.setdefault(only, value)

    rerooted = Hierarchy(pairs)
    rerooted.support = {node: value for node, value in support.items() if rerooted.children[node]}
    rerooted.lengths = lengths
    logging.info(f"Rooted the tree on the branch above {', '.join(sorted(outgroup))}")
    return rerooted
