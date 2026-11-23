"""Collapsing a genome tree into a scheme.

Every bifurcation of a genome tree is far more nodes than a canSNP scheme has,
and many are poorly supported or rest on a few substitutions. ``collapse``
removes internal nodes whose support, branch length or number of genomes falls
below a threshold, joining each one into its parent, and places every genome at
the deepest node that remains. The result is a hierarchy and a genome map that
the extract subcommand reads directly.
"""

import logging
import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

from cansnpscorer.hierarchy import Hierarchy
from cansnpscorer.trees.rooting import outgroup_leaves, reroot


@dataclass
class Collapsed:
    """A genome tree reduced to the nodes of a scheme."""

    hierarchy: Hierarchy
    genomes: dict[str, str]
    """Each genome and the node it sits at."""
    source: dict[str, str]
    """Each node and the node of the genome tree it was, before any renaming."""
    removed: dict[str, str] = field(default_factory=dict)
    """Each internal node of the tree that was removed, and why."""

    def genome_count(self) -> dict[str, int]:
        """Genomes at or below each node."""
        own: dict[str, int] = {}
        for node in self.genomes.values():
            own[node] = own.get(node, 0) + 1
        count: dict[str, int] = {}
        for node in reversed(self.hierarchy.nodes):
            count[node] = own.get(node, 0) + sum(count[child] for child in self.hierarchy.children[node])
        return count


def collapse(
    tree: Hierarchy,
    min_support: float | None = None,
    min_genomes: int = 2,
    min_length: float | None = None,
    members: Mapping[str, Sequence[str]] | None = None,
    outgroups: Sequence[str] = (),
) -> Collapsed:
    """Remove weak internal nodes of a genome tree and place its genomes.

    An internal node other than the root is kept when its support is at
    least ``min_support``, the branch above it is at least ``min_length``,
    and at least ``min_genomes`` genomes lie below it. A node without a
    support value or branch length is not judged on it. Every leaf is a
    genome and sits at the deepest kept node above it; the genomes in
    ``members`` sit wherever the leaf they were folded into sits.

    A kept node that would hold exactly the genomes of its parent, because it
    is the parent's only child and the parent holds no genome of its own, is
    joined into the parent too: no SNP could tell the two apart.

    The leaves named in ``outgroups`` are placed at no node. They stay in the
    alignment as genomes outside every node, which is what gives the root and
    its children a polarity. The tree is first rooted on the branch above the
    outgroups, as ``reroot`` does, so the root of the result is the ingroup's.
    """
    members = members or {}
    outgroup = set(outgroup_leaves(tree, outgroups))
    tree = reroot(tree, sorted(outgroup))
    nodes = tree.nodes
    root = tree.root
    weight = {leaf: 0 if leaf in outgroup else 1 + len(members.get(leaf, [])) for leaf in tree.leaves}
    below: dict[str, int] = {}
    for node in reversed(nodes):
        below[node] = weight.get(node, 0) + sum(below[child] for child in tree.children[node])

    removed: dict[str, str] = {}
    kept: set[str] = {root}
    for node in nodes:
        if node == root or not tree.children[node]:
            continue
        reasons = []
        if min_support is not None and node in tree.support and tree.support[node] < min_support:
            reasons.append(f"support {tree.support[node]:g}")
        if min_length is not None and node in tree.lengths and tree.lengths[node] < min_length:
            reasons.append(f"branch length {tree.lengths[node]:g}")
        if below[node] < min_genomes:
            reasons.append(f"{below[node]} genomes")
        if reasons:
            removed[node] = ", ".join(reasons)
        else:
            kept.add(node)

    # The deepest kept node above each node, found from the root down.
    kept_above: dict[str, str] = {}
    for node in nodes[1:]:
        parent = tree.parents[node]
        # Every node after the first, the root, has a parent.
        assert parent is not None
        kept_above[node] = parent if parent in kept else kept_above[parent]

    parent_of = {node: (None if node == root else kept_above[node]) for node in nodes if node in kept}
    genomes: dict[str, str] = {}
    for leaf in tree.leaves:
        if leaf == root or leaf in outgroup:
            continue
        node = kept_above[leaf]
        genomes[leaf] = node
        for member in members.get(leaf, []):
            genomes[member] = node

    # Join a node that holds exactly its parent's genomes into the parent:
    # its parent has no other child and holds no genome of its own. Joining
    # one node does not change whether another is joined, so the nodes are
    # judged on the tree of kept nodes as it stands, from the root down, and
    # each is joined into the highest node of its run.
    child_count: dict[str, int] = {}
    for parent in parent_of.values():
        if parent is not None:
            child_count[parent] = child_count.get(parent, 0) + 1
    holding = set(genomes.values())
    joined_into: dict[str, str] = {}
    for node, parent in parent_of.items():
        if parent is not None and child_count[parent] == 1 and parent not in holding:
            joined_into[node] = joined_into.get(parent, parent)
            removed[node] = "holds the same genomes as its parent"
    if joined_into:
        parent_of = {
            node: (None if parent is None else joined_into.get(parent, parent))
            for node, parent in parent_of.items()
            if node not in joined_into
        }
        genomes = {genome: joined_into.get(node, node) for genome, node in genomes.items()}

    hierarchy = Hierarchy(list(parent_of.items()))
    logging.info(
        f"Kept {len(hierarchy)} of {len(tree) - len(tree.leaves)} internal nodes, and placed "
        f"{len(genomes)} genomes"
    )
    return Collapsed(
        hierarchy=hierarchy, genomes=genomes, source={node: node for node in hierarchy.nodes}, removed=removed
    )


def rename(collapsed: Collapsed, prefix: str) -> Collapsed:
    """Give the nodes hierarchical names: the root ``prefix``, its children
    ``prefix.1``, ``prefix.2`` and so on, largest first.

    The names of a genome tree are hashes or generated numbers, which say
    nothing and change when genomes are added. Children with as many genomes
    are ordered by the first of their genome names in sort order, so that
    the names depend on the genomes and not on the order the tree file lists
    its branches in or the names it gives its nodes.
    """
    if not prefix or ";" in prefix or "\t" in prefix:
        raise ValueError("the prefix must be a non-empty name without semicolons or tabs")
    count = collapsed.genome_count()
    tree = collapsed.hierarchy
    first: dict[str, str] = {}
    for genome, node in collapsed.genomes.items():
        if node not in first or genome < first[node]:
            first[node] = genome
    for node in reversed(tree.nodes):
        for child in tree.children[node]:
            if child in first and (node not in first or first[child] < first[node]):
                first[node] = first[child]
    names = {tree.root: prefix}
    for node in tree.nodes:
        # A node without genomes, which only a min_genomes below 1 keeps, last.
        ordered = sorted(
            tree.children[node],
            key=lambda child: (-count[child], child not in first, first.get(child, ""), child),
        )
        for index, child in enumerate(ordered, 1):
            names[child] = f"{names[node]}.{index}"

    hierarchy = Hierarchy(
        [(names[node], None if parent is None else names[parent]) for node, parent in tree.parents.items()]
    )
    return Collapsed(
        hierarchy=hierarchy,
        genomes={genome: names[node] for genome, node in collapsed.genomes.items()},
        source={names[node]: collapsed.source[node] for node in tree.nodes},
        removed=collapsed.removed,
    )


def write_collapsed(collapsed: Collapsed, tree: Hierarchy, outdir: str, prefix: str) -> list[str]:
    """Write the hierarchy, the genome map and a node table."""
    os.makedirs(outdir, exist_ok=True)
    hierarchy = collapsed.hierarchy
    count = collapsed.genome_count()
    paths = [
        os.path.join(outdir, f"{prefix}_{name}") for name in ("hierarchy.tsv", "genomes.tsv", "nodes.tsv")
    ]

    with open(paths[0], "w", newline="\n") as handle:
        handle.write("child\tparent\n")
        handle.writelines(f"{node}\t{hierarchy.parents[node] or ''}\n" for node in hierarchy.nodes)

    with open(paths[1], "w", newline="\n") as handle:
        handle.writelines(
            f"{genome}\t{node}\n"
            for genome, node in sorted(collapsed.genomes.items(), key=lambda item: (item[1], item[0]))
        )

    with open(paths[2], "w", newline="\n") as handle:
        handle.write("node\tparent\tgenomes\tsupport\tbranch_length\ttree_node\n")
        for node in hierarchy.nodes:
            source = collapsed.source[node]
            support = tree.support.get(source)
            length = tree.lengths.get(source)
            handle.write(
                "\t".join(
                    [
                        node,
                        hierarchy.parents[node] or "",
                        str(count[node]),
                        "NA" if support is None else f"{support:g}",
                        "NA" if length is None else f"{length:g}",
                        source,
                    ]
                )
                + "\n"
            )
    return paths
