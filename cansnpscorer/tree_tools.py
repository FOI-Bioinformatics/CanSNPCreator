"""Genome trees: collapsing them into a scheme, and placing genomes on one.

A genome tree, as RepGenR2 builds it, has one leaf per genome and one internal
node per split the tree builder found. Two things follow from it.

**Collapsing.** Every bifurcation of a genome tree is far more nodes than a
canSNP scheme has, and many are poorly supported or rest on a few substitutions.
``collapse`` removes internal nodes whose support, branch length or number of
genomes falls below a threshold, joining each one into its parent, and places
every genome at the deepest node that remains. The result is a hierarchy and a
genome map that the extract subcommand reads directly.

**Placing.** To check an existing scheme against genomes, each genome must sit
at a node of that scheme. Typing the genomes with the scheme itself places them
by the very SNPs under test, so a genome lacking the derived base of its own
branch is typed elsewhere and never seen to lack it. ``place`` uses each
scheme SNP only to identify which clade of the tree the scheme's node means:
the clade whose genomes best match the genomes carrying the derived base. The
genomes are then placed by the topology of the tree, so a genome that sits
inside the clade but lacks the derived base stays inside, where the audit
reports it.

The match is by Jaccard index between the genomes carrying the derived base and
the genomes on one side of a branch of the tree, counting only genomes that
have a base at the site. Both sides of every branch are candidates, so the
match does not depend on where the tree is rooted: a midpoint root that falls
inside a clade of the scheme still leaves that clade on one side of a branch.
A clade of two genomes in which one lacks the derived base matches its other
genome equally well, so placement is robust for nodes with several genomes and
weak for nodes with one or two; the report gives the index of every match.
"""

import csv
import logging
import os
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, List, Mapping, Optional, Sequence, Set, Tuple

from cansnpscorer.findings import CHECK, MODIFY, Finding, problem, recommendation
from cansnpscorer.helpers.file_utils import read_hc_file, recognize_hc_file_format
from cansnpscorer.helpers.names import resolve_names
from cansnpscorer.hierarchy import Hierarchy, HierarchyError
from cansnpscorer.model import AlignmentMap, Reference
from cansnpscorer.scheme import Scheme

BASES = frozenset("ACGT")
NODE_TABLE_COLUMNS = ("node", "support", "branch_length")


# -- reading a genome tree ---------------------------------------------------


def read_genome_tree(path: str, nodes_file: Optional[str] = None) -> Hierarchy:
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
    tree = read_hc_file(path, recognize_hc_file_format(path))
    if nodes_file:
        with open(nodes_file) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = [name for name in NODE_TABLE_COLUMNS
                       if name not in (reader.fieldnames or [])]
            if missing:
                raise ValueError(
                    f"{nodes_file}: lacks the columns {', '.join(missing)}; it should be "
                    f"RepGenR2's tree2tax_nodes.tsv")
            for row in reader:
                node = row["node"].strip()
                if node not in tree:
                    continue
                for column, target in (("support", tree.support),
                                       ("branch_length", tree.lengths)):
                    value = (row.get(column) or "").strip()
                    if value and value.upper() != "NA":
                        target[node] = float(value)
    logging.info(
        f"Read a genome tree of {len(tree.leaves)} genomes and "
        f"{len(tree) - len(tree.leaves)} internal nodes, {len(tree.support)} with a support "
        f"value and {len(tree.lengths)} with a branch length")
    return tree


def read_members(path: str, tree: Hierarchy) -> Dict[str, List[str]]:
    """Genomes a leaf of the tree stands for, from RepGenR2's genomes_map.tsv.

    That file maps an accession to a leaf. The leaf's own accession is one of
    its rows and is left out here; the other rows are the genomes dereplication
    folded into that leaf.
    """
    leaves = set(tree.leaves)
    members: Dict[str, List[str]] = {}
    with open(path) as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if not line.strip() or line.startswith("#"):
                continue
            if len(fields) != 2:
                raise ValueError(f"{path}: line {number} should hold an accession and a leaf")
            genome, leaf = fields[0].strip(), fields[1].strip()
            if leaf not in leaves:
                continue
            own, _unmatched, _ambiguous = resolve_names([genome], [leaf])
            if own:
                continue
            members.setdefault(leaf, []).append(genome)
    logging.info(
        f"Read {sum(len(names) for names in members.values())} genomes folded into "
        f"{len(members)} leaves of the tree")
    return members


# -- collapsing ----------------------------------------------------------------


@dataclass
class Collapsed:
    """A genome tree reduced to the nodes of a scheme."""

    hierarchy: Hierarchy
    genomes: Dict[str, str]
    """Each genome and the node it sits at."""
    source: Dict[str, str]
    """Each node and the node of the genome tree it was, before any renaming."""
    removed: Dict[str, str] = field(default_factory=dict)
    """Each internal node of the tree that was removed, and why."""

    def genome_count(self) -> Dict[str, int]:
        """Genomes at or below each node."""
        own: Dict[str, int] = {}
        for node in self.genomes.values():
            own[node] = own.get(node, 0) + 1
        count: Dict[str, int] = {}
        for node in reversed(self.hierarchy.nodes):
            count[node] = own.get(node, 0) + sum(
                count[child] for child in self.hierarchy.children[node])
        return count


def collapse(tree: Hierarchy, min_support: Optional[float] = None, min_genomes: int = 2,
             min_length: Optional[float] = None,
             members: Optional[Mapping[str, Sequence[str]]] = None,
             outgroups: Sequence[str] = ()) -> Collapsed:
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
    its children a polarity; the root of the result is then the ingroup's.
    """
    members = members or {}
    unknown = sorted(set(outgroups) - set(tree.leaves))
    if unknown:
        raise HierarchyError(f"the tree has no leaf named {', '.join(unknown)}")
    outgroup = set(outgroups)
    weight = {leaf: 0 if leaf in outgroup else 1 + len(members.get(leaf, []))
              for leaf in tree.leaves}
    below: Dict[str, int] = {}
    for node in reversed(tree.nodes):
        below[node] = weight.get(node, 0) + sum(below[child] for child in tree.children[node])

    removed: Dict[str, str] = {}
    kept: Set[str] = {tree.root}
    for node in tree.nodes:
        if node == tree.root or not tree.children[node]:
            continue
        reasons = []
        if min_support is not None and node in tree.support \
                and tree.support[node] < min_support:
            reasons.append(f"support {tree.support[node]:g}")
        if min_length is not None and node in tree.lengths and tree.lengths[node] < min_length:
            reasons.append(f"branch length {tree.lengths[node]:g}")
        if below[node] < min_genomes:
            reasons.append(f"{below[node]} genomes")
        if reasons:
            removed[node] = ", ".join(reasons)
        else:
            kept.add(node)

    def kept_above(node: str) -> str:
        parent = tree.parents[node]
        while parent not in kept:
            parent = tree.parents[parent]
        return parent

    parent_of = {node: (None if node == tree.root else kept_above(node))
                 for node in tree.nodes if node in kept}
    genomes: Dict[str, str] = {}
    for leaf in tree.leaves:
        if leaf == tree.root or leaf in outgroup:
            continue
        node = kept_above(leaf)
        genomes[leaf] = node
        for member in members.get(leaf, []):
            genomes[member] = node

    # Join a node that holds exactly its parent's genomes into the parent.
    while True:
        children: Dict[str, List[str]] = {}
        for node, parent in parent_of.items():
            if parent is not None:
                children.setdefault(parent, []).append(node)
        holding = set(genomes.values())
        joined = next(
            (node for node, parent in parent_of.items()
             if parent is not None and len(children.get(parent, [])) == 1
             and parent not in holding),
            None)
        if joined is None:
            break
        parent = parent_of.pop(joined)
        for node, above in parent_of.items():
            if above == joined:
                parent_of[node] = parent
        for genome, node in genomes.items():
            if node == joined:
                genomes[genome] = parent
        removed[joined] = "holds the same genomes as its parent"

    hierarchy = Hierarchy(list(parent_of.items()))
    logging.info(
        f"Kept {len(hierarchy)} of {len(tree) - len(tree.leaves)} internal nodes, and placed "
        f"{len(genomes)} genomes")
    return Collapsed(hierarchy=hierarchy, genomes=genomes,
                     source={node: node for node in hierarchy.nodes}, removed=removed)


def rename(collapsed: Collapsed, prefix: str) -> Collapsed:
    """Give the nodes hierarchical names: the root ``prefix``, its children
    ``prefix.1``, ``prefix.2`` and so on, largest first.

    The names of a genome tree are hashes or generated numbers, which say
    nothing and change when genomes are added.
    """
    if not prefix or ";" in prefix or "\t" in prefix:
        raise ValueError("the prefix must be a non-empty name without semicolons or tabs")
    count = collapsed.genome_count()
    tree = collapsed.hierarchy
    names = {tree.root: prefix}
    for node in tree.nodes:
        ordered = sorted(tree.children[node], key=lambda child: (-count[child], child))
        for index, child in enumerate(ordered, 1):
            names[child] = f"{names[node]}.{index}"

    hierarchy = Hierarchy([(names[node], None if parent is None else names[parent])
                           for node, parent in tree.parents.items()])
    return Collapsed(
        hierarchy=hierarchy,
        genomes={genome: names[node] for genome, node in collapsed.genomes.items()},
        source={names[node]: collapsed.source[node] for node in tree.nodes},
        removed=collapsed.removed)


def write_collapsed(collapsed: Collapsed, tree: Hierarchy, outdir: str,
                    prefix: str) -> List[str]:
    """Write the hierarchy, the genome map and a node table."""
    os.makedirs(outdir, exist_ok=True)
    hierarchy = collapsed.hierarchy
    count = collapsed.genome_count()
    paths = [os.path.join(outdir, f"{prefix}_{name}")
             for name in ("hierarchy.tsv", "genomes.tsv", "nodes.tsv")]

    with open(paths[0], "w", newline="\n") as handle:
        handle.write("child\tparent\n")
        for node in hierarchy.nodes:
            handle.write(f"{node}\t{hierarchy.parents[node] or ''}\n")

    with open(paths[1], "w", newline="\n") as handle:
        for genome, node in sorted(collapsed.genomes.items(), key=lambda item: (item[1], item[0])):
            handle.write(f"{genome}\t{node}\n")

    with open(paths[2], "w", newline="\n") as handle:
        handle.write("node\tparent\tgenomes\tsupport\tbranch_length\ttree_node\n")
        for node in hierarchy.nodes:
            source = collapsed.source[node]
            support = tree.support.get(source)
            length = tree.lengths.get(source)
            handle.write("\t".join([
                node, hierarchy.parents[node] or "", str(count[node]),
                "NA" if support is None else f"{support:g}",
                "NA" if length is None else f"{length:g}", source]) + "\n")
    return paths


# -- placing genomes on an existing scheme -----------------------------------


@dataclass
class Match:
    """The clade of the genome tree a node of the scheme was matched to."""

    node: str
    clade: Optional[str] = None
    jaccard: float = 0.0
    carriers: int = 0
    shared: int = 0
    clade_only: int = 0
    carriers_only: int = 0
    genomes: FrozenSet[str] = frozenset()
    reason: str = ""
    outside: bool = False
    """Whether the match is the genomes outside ``clade`` rather than inside it."""
    group: FrozenSet[str] = frozenset()
    """The genomes the node was matched to."""

    @property
    def described(self) -> str:
        return f"the genomes outside {self.clade}" if self.outside else str(self.clade)

    HEADER = ("node", "tree_node", "side", "jaccard", "carriers", "shared", "clade_only",
              "carriers_only", "genomes_placed", "note")


@dataclass
class Placement:
    """Where the genomes sit in an existing scheme, read from a genome tree."""

    clusters: Dict[str, List[str]]
    """The genomes placed at each node of the scheme, by alignment record name."""
    matches: List[Match]
    findings: List[Finding]


def _clades(tree: Hierarchy, records: Mapping[str, str]) -> Dict[str, FrozenSet[str]]:
    """The alignment records at or below each node of the tree."""
    clades: Dict[str, FrozenSet[str]] = {}
    for node in reversed(tree.nodes):
        own = {records[node]} if node in records else set()
        for child in tree.children[node]:
            own |= clades[child]
        clades[node] = frozenset(own)
    return clades


def place(scheme: Scheme, tree: Hierarchy, alignment: Mapping[str, str], reference: Reference,
          genome: Optional[str] = None, min_match: float = 0.5,
          members: Optional[Mapping[str, Sequence[str]]] = None) -> Placement:
    """Place the genomes of a tree on the nodes of an existing scheme.

    Args:
        scheme: The scheme, as read from its source files.
        tree: A rooted genome tree whose leaves are the alignment's genomes.
        alignment: The whole genome alignment, reference first.
        reference: The reference genome the scheme and the alignment share.
        genome: Which of the scheme's reference genomes that is, when it names
            more than one.
        min_match: The lowest Jaccard index accepted as a match.
        members: Genomes folded into a leaf of the tree, which are placed with
            it when the alignment holds them.
    """
    findings: List[Finding] = []
    if scheme.hierarchy is None:
        return Placement({}, [], [problem(
            CHECK, "placement", "", "the tree file of the scheme is not a single rooted tree")])

    mapping = AlignmentMap(next(iter(alignment.values())), reference)
    records = list(alignment)
    matched, unmatched, ambiguous = resolve_names(tree.leaves, records)
    if unmatched or ambiguous:
        findings.append(recommendation(
            CHECK, "tree-genomes", "",
            f"{len(unmatched) + len(ambiguous)} leaves of the tree match no single genome of "
            f"the alignment and are left out: {', '.join(sorted(unmatched + ambiguous)[:5])}"))
    leaf_record = dict(matched)

    member_record: Dict[str, str] = {}
    for leaf, names in (members or {}).items():
        if leaf not in leaf_record:
            continue
        found, _missing, _ambiguous = resolve_names(names, records)
        for record in found.values():
            member_record[record] = leaf

    clades = _clades(tree, leaf_record)
    in_tree = clades[tree.root]

    matches: Dict[str, Match] = {}
    snps = {}
    for snp in scheme.snps:
        if genome is not None and snp.genome != genome:
            continue
        snps.setdefault(snp.node, snp)

    for node in scheme.hierarchy.nodes:
        snp = snps.get(node)
        if snp is None:
            matches[node] = Match(node, reason="no SNP for this node")
            continue
        if not 1 <= snp.position <= reference.length:
            matches[node] = Match(node, reason="the position is outside the reference")
            continue
        column = mapping.to_column(snp.position - 1)
        if column is None:
            matches[node] = Match(node, reason="the position is not aligned")
            continue

        called: Set[str] = set()
        carriers: Set[str] = set()
        for record in in_tree:
            base = alignment[record][column].upper()
            if base in BASES:
                called.add(record)
                if base == snp.derived.upper():
                    carriers.add(record)
        if not carriers:
            matches[node] = Match(node, reason="no genome of the tree carries the derived base")
            continue

        # Count, for every node of the tree, the genomes below it that have a
        # base here and those that carry the derived one. Each branch splits
        # the genomes into the clade below it and the rest; both sides are
        # candidates, so the match does not depend on where the tree is rooted.
        n_called = len(called)
        n_carriers = len(carriers)
        below_called: Dict[str, int] = {}
        below_carriers: Dict[str, int] = {}
        for tree_node in reversed(tree.nodes):
            record = leaf_record.get(tree_node)
            own_called = 1 if record in called else 0
            own_carrier = 1 if record in carriers else 0
            below_called[tree_node] = own_called + sum(
                below_called[child] for child in tree.children[tree_node])
            below_carriers[tree_node] = own_carrier + sum(
                below_carriers[child] for child in tree.children[tree_node])

        best: Optional[Tuple[float, int, int, str]] = None
        for tree_node in tree.nodes:
            inside_called, inside_carriers = below_called[tree_node], below_carriers[tree_node]
            sides = [(inside_called, inside_carriers, 0)]
            if tree_node != tree.root:
                sides.append((n_called - inside_called, n_carriers - inside_carriers, 1))
            for side_called, shared, outside in sides:
                if not shared:
                    continue
                index = shared / (side_called + n_carriers - shared)
                size = (len(clades[tree_node]) if not outside
                        else len(in_tree) - len(clades[tree_node]))
                # The higher index first; then the side that is a clade of the
                # tree as rooted; then the smaller side.
                key = (-index, outside, size, tree_node)
                if best is None or key < best:
                    best = key

        index, outside, _size, clade = -best[0], bool(best[1]), best[2], best[3]
        group = in_tree - clades[clade] if outside else clades[clade]
        genomes = group & called
        match = Match(node, clade, index, len(carriers), len(genomes & carriers),
                      len(genomes - carriers), len(carriers - genomes),
                      outside=outside, group=group)
        if index < min_match:
            match.reason = f"the best clade matches with a Jaccard index of {index:.2f} only"
            match.clade = None
        elif outside:
            match.reason = (f"matches the genomes outside the tree's clade {clade}, so the tree "
                            f"is rooted inside this node's clade")
        matches[node] = match

    for node, match in matches.items():
        if match.clade is None:
            findings.append(recommendation(
                CHECK, "unplaced-node", node,
                f"the node could not be matched to a clade of the tree: {match.reason}. No "
                f"genome is placed at it"))
        elif match.jaccard < 1.0:
            findings.append(recommendation(
                CHECK, "partial-match", node,
                f"the node matches {match.described} with a Jaccard index of "
                f"{match.jaccard:.2f}: {match.clade_only} genomes of the clade lack the derived "
                f"base, and {match.carriers_only} genomes outside it carry it"))

    rerooted = sorted(node for node, match in matches.items() if match.clade and match.outside)
    if rerooted:
        findings.append(recommendation(
            CHECK, "tree-root", ", ".join(rerooted[:5]),
            f"{len(rerooted)} nodes match the genomes outside a clade of the tree rather than a "
            f"clade, so the tree is rooted inside them. The placement does not depend on it, "
            f"but a tree rooted on an outgroup avoids the question"))

    for node in scheme.hierarchy.nodes:
        parent = scheme.hierarchy.parents[node]
        child_match, parent_match = matches[node], matches.get(parent) if parent else None
        if child_match.clade and parent_match and parent_match.clade \
                and not child_match.group <= parent_match.group:
            findings.append(problem(
                MODIFY, "placement-order", node,
                f"the node matches {child_match.described}, which is not inside "
                f"{parent_match.described}, which its parent {parent} matches. The tree and the "
                f"scheme disagree on where this branch belongs"))

    depth = {node: len(scheme.hierarchy.path(node)) for node in scheme.hierarchy.nodes}
    matched_nodes = [node for node, match in matches.items() if match.clade]
    clusters: Dict[str, List[str]] = {}
    split: List[str] = []
    placed_at: Dict[str, str] = {}
    for record in sorted(in_tree):
        holding = [node for node in matched_nodes if record in matches[node].group]
        if not holding:
            continue
        deepest = max(holding, key=lambda node: (depth[node], matches[node].jaccard))
        path = set(scheme.hierarchy.path(deepest))
        if any(node not in path for node in holding):
            split.append(record)
        placed_at[record] = deepest
    for record, leaf in member_record.items():
        if record not in placed_at and leaf_record[leaf] in placed_at:
            placed_at[record] = placed_at[leaf_record[leaf]]
    for record, node in placed_at.items():
        clusters.setdefault(node, []).append(record)

    if split:
        findings.append(problem(
            CHECK, "placement-conflict", "",
            f"{len(split)} genomes lie in the clades of nodes on different branches of the "
            f"scheme, and were placed at the deepest of them: {', '.join(split[:5])}"))

    for match in matches.values():
        match.genomes = frozenset(clusters.get(match.node, []))

    logging.info(
        f"Matched {len(matched_nodes)} of {len(matches)} nodes of the scheme to clades of the "
        f"tree, and placed {len(placed_at)} of {len(in_tree)} genomes")
    return Placement(clusters, list(matches.values()), findings)


def write_placement(placement: Placement, genome_map: str, report: Optional[str] = None) -> None:
    """Write the genome map, and a report of how each node was matched."""
    with open(genome_map, "w", newline="\n") as handle:
        for node, genomes in sorted(placement.clusters.items()):
            for genome in sorted(genomes):
                handle.write(f"{genome}\t{node}\n")
    if report:
        with open(report, "w", newline="\n") as handle:
            handle.write("\t".join(Match.HEADER) + "\n")
            for match in placement.matches:
                handle.write("\t".join([
                    match.node, match.clade or "",
                    "" if not match.clade else ("outside" if match.outside else "inside"),
                    f"{match.jaccard:.3f}", str(match.carriers),
                    str(match.shared), str(match.clade_only), str(match.carriers_only),
                    str(len(match.genomes)), match.reason]) + "\n")

