"""Reading a genome tree and the genomes of its clusters."""

import csv
import logging

from cansnpscorer.helpers.file_utils import hierarchy_format, open_text, read_hierarchy
from cansnpscorer.helpers.names import resolve_names
from cansnpscorer.hierarchy import Hierarchy

NODE_TABLE_COLUMNS = ("node", "support", "branch_length")


def read_genome_tree(path: str, nodes_file: str | None = None) -> Hierarchy:
    """Read a genome tree whose leaves are genomes.

    Args:
        path: A Newick tree, or a child and parent table such as RepGenR2's
            tree2tax.tsv.
        nodes_file: RepGenR2's tree2tax_nodes.tsv, whose support and
            branch_length columns are added to the tree. A Newick tree carries
            both itself.

    Raises:
        HierarchyError: If the tree cannot be read.
        ValueError: If the node table lacks its columns.
    """
    if hierarchy_format(path) == "newick":
        # A genome tree's numeric internal labels are support values.
        with open_text(path) as handle:
            tree = Hierarchy.from_newick(handle.read(), support_labels=True)
    else:
        tree = read_hierarchy(path)
    if nodes_file:
        with open_text(nodes_file) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = [name for name in NODE_TABLE_COLUMNS if name not in (reader.fieldnames or [])]
            if missing:
                raise ValueError(
                    f"{nodes_file}: lacks the columns {', '.join(missing)}; it should be "
                    f"RepGenR2's tree2tax_nodes.tsv"
                )
            for row in reader:
                node = row["node"].strip()
                if node not in tree:
                    continue
                for column, target in (("support", tree.support), ("branch_length", tree.lengths)):
                    value = (row.get(column) or "").strip()
                    if value and value.upper() != "NA":
                        target[node] = float(value)
        # As for Newick, supports written as percentages are read as fractions.
        if any(value > 1.0 for value in tree.support.values()):
            tree.support = {node: value / 100.0 for node, value in tree.support.items()}
    logging.info(
        f"Read a genome tree of {len(tree.leaves)} genomes and "
        f"{len(tree) - len(tree.leaves)} internal nodes, {len(tree.support)} with a support "
        f"value and {len(tree.lengths)} with a branch length"
    )
    return tree


def read_members(path: str, tree: Hierarchy) -> dict[str, list[str]]:
    """Genomes a leaf of the tree stands for, from RepGenR2's genomes_map.tsv.

    That file maps an accession to a leaf. The leaf's own accession is one of
    its rows and is left out here; the other rows are the genomes dereplication
    folded into that leaf.
    """
    leaves = set(tree.leaves)
    members: dict[str, list[str]] = {}
    foreign: set[str] = set()
    placed: dict[str, str] = {}
    twice: list[str] = []
    rows: list[tuple[str, str]] = []
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if not line.strip() or line.startswith("#"):
                continue
            if len(fields) != 2:
                raise ValueError(f"{path}: line {number} should hold an accession and a leaf")
            genome, leaf = fields[0].strip(), fields[1].strip()
            if leaf not in leaves:
                foreign.add(leaf)
                continue
            rows.append((genome, leaf))
    # Matched against the leaves once, not row by row.
    is_leaf, _unmatched, _ambiguous = resolve_names(sorted({genome for genome, _ in rows}), tree.leaves)
    for genome, leaf in rows:
        own, _unmatched, _ambiguous = resolve_names([genome], [leaf])
        if own:
            continue
        # A genome that is a leaf itself, as in a tree of all genomes, sits
        # at its own leaf; one listed under two leaves keeps the first.
        if genome in is_leaf:
            continue
        if genome in placed:
            if placed[genome] != leaf:
                twice.append(genome)
            continue
        placed[genome] = leaf
        members.setdefault(leaf, []).append(genome)
    if twice:
        logging.warning(
            f"{path} places {len(twice)} genomes under two leaves; each was kept under the "
            f"first: {', '.join(sorted(twice)[:3])}"
        )
    if foreign:
        logging.warning(
            f"{path} places genomes at {len(foreign)} leaves the tree does not have, and they "
            f"were left out: {', '.join(sorted(foreign)[:3])}. The file may come from another "
            f"tree; RepGenR2 writes it with tree2tax, which must be rerun after phylo"
        )
    logging.info(
        f"Read {sum(len(names) for names in members.values())} genomes folded into "
        f"{len(members)} leaves of the tree"
    )
    return members
