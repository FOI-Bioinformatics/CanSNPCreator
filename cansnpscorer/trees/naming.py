"""Naming the nodes of a hierarchy after the labelled genomes they hold."""

import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from cansnpscorer.hierarchy import Hierarchy

#: Characters no node name may hold, as CanSNPer2-data requires.
NAME_REFUSED = (";", "\t", "\n", "\r")


@dataclass
class Named:
    """A hierarchy whose nodes are named after the labelled genomes they hold."""

    hierarchy: Hierarchy
    genomes: dict[str, str]  # each genome -> its node, under the new names
    names: dict[str, str]  # each node's old name -> its new name
    named_by: dict[str, str]  # each node's old name -> root, label or position
    labels: dict[str, str]  # each node named after a label -> that label
    below: dict[str, int]  # each node's old name -> the genomes at or below it
    labelled: dict[str, int]  # each node's old name -> the labelled genomes among them
    unplaced: list[str]  # labels that no node holds apart from another label
    split: dict[str, list[str]]  # labels named at more than one node -> those nodes
    missing: list[str]  # labelled genomes the genome map does not place


def name_from_labels(
    hierarchy: Hierarchy, clusters: Mapping[str, Sequence[str]], labels: Mapping[str, str]
) -> Named:
    """Name the nodes of a hierarchy after the labels of the genomes they hold.

    A node takes a label when every labelled genome at or below it carries that
    label, and of the nodes that do, the shallowest takes it, so the name covers
    as much of the clade as the labels allow. The root keeps its name, and a
    named node becomes ``<root>.<label>``. Every other node keeps its last name
    segment under its parent's new name, so the nodes below a named node extend
    its name: with ``Bsu.1.1.2`` named ``Bsu.bv2``, ``Bsu.1.1.2.1`` becomes
    ``Bsu.bv2.1``.

    Most genomes carry no label, and need none. A label that shares every node
    with another label is not placed, and a label that ends up at two nodes,
    because its genomes do not form a clade, names both: the node with more
    genomes takes the plain name and the other a numbered one. Both cases are
    reported, since they are findings about the labels rather than errors.
    """
    for genome, label in labels.items():
        if not label or any(character in label for character in NAME_REFUSED):
            raise ValueError(
                f"the label {label!r} of {genome} cannot name a node: it must be non-empty and "
                f"hold no semicolon, tab or line break"
            )

    placed = {genome for genomes in clusters.values() for genome in genomes}
    missing = sorted(set(labels) - placed)
    labelled = {genome: label for genome, label in labels.items() if genome in placed}

    root = hierarchy.root
    isolates = hierarchy.isolates(clusters)
    below = {node: len(isolates[node]) for node in hierarchy.nodes}
    held = {node: {labelled[g] for g in isolates[node] if g in labelled} for node in hierarchy.nodes}

    # The shallowest node holding one label alone: its parent holds another
    # label too, or is the root, which keeps its own name.
    label_at = {}
    for node in hierarchy.nodes:
        if node == root or len(held[node]) != 1:
            continue
        parent = hierarchy.parents[node]
        # Only the root has no parent, and it is skipped above.
        assert parent is not None
        if parent == root or held[parent] != held[node]:
            label_at[node] = next(iter(held[node]))

    used = {root}

    def unique(base: str) -> str:
        name, number = base, 2
        while name in used:
            name, number = f"{base}-{number}", number + 1
        used.add(name)
        return name

    # Labelled names first, the larger clade taking the plain name, so that
    # which node of a split label gets it does not depend on traversal order.
    by_label: dict[str, list[str]] = {}
    for node, label in label_at.items():
        by_label.setdefault(label, []).append(node)
    reserved = {}
    for label in sorted(by_label):
        for node in sorted(by_label[label], key=lambda name: (-below[name], name)):
            reserved[node] = unique(f"{root}.{label}")

    names = {root: root}
    named_by = {root: "root"}
    for node in hierarchy.nodes:
        if node == root:
            continue
        if node in reserved:
            names[node], named_by[node] = reserved[node], "label"
            continue
        # Only a name that extends its parent's, as collapse --rename writes them,
        # follows the parent to its new name: Bsu.1.1.2.1 under a renamed
        # Bsu.1.1.2 becomes Bsu.bv2.1. A name that is its own, such as S.5 under
        # B.2, is kept as it is.
        old_parent = hierarchy.parents[node]
        # Only the root has no parent, and it is skipped above.
        assert old_parent is not None
        if node.startswith(old_parent + ".") and names[old_parent] != old_parent:
            names[node] = unique(names[old_parent] + node[len(old_parent) :])
        else:
            used.add(node)
            names[node] = node
        named_by[node] = "position"

    renamed = Hierarchy(
        [
            (names[node], None if parent is None else names[parent])
            for node, parent in hierarchy.parents.items()
        ]
    )
    return Named(
        hierarchy=renamed,
        genomes={genome: names[node] for node, genomes in clusters.items() for genome in genomes},
        names=names,
        named_by=named_by,
        labels=label_at,
        below=below,
        labelled={node: sum(1 for g in isolates[node] if g in labelled) for node in hierarchy.nodes},
        unplaced=sorted({*labelled.values()} - set(by_label)),
        split={label: sorted(nodes) for label, nodes in by_label.items() if len(nodes) > 1},
        missing=missing,
    )


def write_named(named: Named, old: Hierarchy, outdir: str, prefix: str) -> list[str]:
    """Write the renamed hierarchy and genome map, and how each node got its name."""
    os.makedirs(outdir, exist_ok=True)
    paths = [
        os.path.join(outdir, f"{prefix}_{name}") for name in ("hierarchy.tsv", "genomes.tsv", "names.tsv")
    ]
    hierarchy = named.hierarchy
    with open(paths[0], "w", newline="\n") as handle:
        handle.write("child\tparent\n")
        handle.writelines(f"{node}\t{hierarchy.parents[node] or ''}\n" for node in hierarchy.nodes)
    with open(paths[1], "w", newline="\n") as handle:
        handle.writelines(
            f"{genome}\t{node}\n"
            for genome, node in sorted(named.genomes.items(), key=lambda item: (item[1], item[0]))
        )
    with open(paths[2], "w", newline="\n") as handle:
        handle.write("old_name\tnew_name\tnamed_by\tlabel\tgenomes\tlabelled\n")
        handle.writelines(
            f"{node}\t{named.names[node]}\t{named.named_by[node]}\t{named.labels.get(node, '')}"
            f"\t{named.below[node]}\t{named.labelled[node]}\n"
            for node in old.nodes
        )
    return paths
