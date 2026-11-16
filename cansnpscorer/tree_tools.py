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
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

from cansnpscorer.findings import CHECK, MODIFY, Finding, problem, recommendation
from cansnpscorer.helpers.file_utils import hierarchy_format, open_text, read_hierarchy
from cansnpscorer.helpers.names import resolve_names
from cansnpscorer.hierarchy import Hierarchy, HierarchyError
from cansnpscorer.model import AlignmentMap, Reference
from cansnpscorer.scheme import Scheme

BASES = frozenset("ACGT")
NODE_TABLE_COLUMNS = ("node", "support", "branch_length")


# -- reading a genome tree ---------------------------------------------------


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
        # As for Newick, supports written as percentages are read as fractions.
        if any(value > 1.0 for value in tree.support.values()):
            tree.support = {node: value / 100.0 for node, value in tree.support.items()}
    logging.info(
        f"Read a genome tree of {len(tree.leaves)} genomes and "
        f"{len(tree) - len(tree.leaves)} internal nodes, {len(tree.support)} with a support "
        f"value and {len(tree.lengths)} with a branch length")
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
    is_leaf, _unmatched, _ambiguous = resolve_names(sorted({genome for genome, _ in rows}),
                                                    tree.leaves)
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
            f"first: {', '.join(sorted(twice)[:3])}")
    if foreign:
        logging.warning(
            f"{path} places genomes at {len(foreign)} leaves the tree does not have, and they "
            f"were left out: {', '.join(sorted(foreign)[:3])}. The file may come from another "
            f"tree; RepGenR2 writes it with tree2tax, which must be rerun after phylo")
    logging.info(
        f"Read {sum(len(names) for names in members.values())} genomes folded into "
        f"{len(members)} leaves of the tree")
    return members


# -- collapsing ----------------------------------------------------------------


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
            count[node] = own.get(node, 0) + sum(
                count[child] for child in self.hierarchy.children[node])
        return count


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
            f"tree; give the full leaf name")
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
    stripped = Hierarchy([(node, None if node == root else tree.parents[node])
                          for node in tree.nodes if node not in gone])
    stripped.support = {node: value for node, value in tree.support.items()
                        if node not in gone}
    # The branch above the new root led to the removed one.
    stripped.lengths = {node: value for node, value in tree.lengths.items()
                        if node not in gone and node != root}
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
        outgroups_below[node] = sum(outgroups_below[child] for child in children) + (
            node in outgroup)

    n_ingroup = leaves_below[root] - len(outgroup)
    if not n_ingroup:
        raise HierarchyError("every leaf of the tree is an outgroup")

    def splits(node: str) -> bool:
        """Whether the clade below the node is the outgroups or the ingroup."""
        return ((leaves_below[node] == outgroups_below[node] == len(outgroup))
                or (outgroups_below[node] == 0 and leaves_below[node] == n_ingroup))

    # The branch that splits the tree into outgroups and ingroup is held by
    # the node below it, on whichever side of the split that node is.
    branch = next((node for node in nodes if node != root and splits(node)), None)
    if branch is None:
        raise HierarchyError(
            f"the outgroups ({', '.join(sorted(outgroup))}) do not form a clade of the tree, "
            f"however it is rooted, so no root separates them from the other genomes")
    if tree.parents[branch] == root and (outgroups_below[branch] == len(outgroup)
                                         or len(tree.children[root]) == 2):
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
    children: dict[str, list[str]] = {}
    for child, parent in pairs:
        if parent is not None:
            children.setdefault(parent, []).append(child)
    old = root
    if old != upper and len(children.get(old, [])) == 1:
        (only,) = children[old]
        parent = next(parent for child, parent in pairs if child == old)
        pairs = [(child, parent if child == only else above) for child, above in pairs
                 if child != old]
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


def collapse(tree: Hierarchy, min_support: float | None = None, min_genomes: int = 2,
             min_length: float | None = None,
             members: Mapping[str, Sequence[str]] | None = None,
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
    its children a polarity. The tree is first rooted on the branch above the
    outgroups, as ``reroot`` does, so the root of the result is the ingroup's.
    """
    members = members or {}
    outgroup = set(outgroup_leaves(tree, outgroups))
    tree = reroot(tree, sorted(outgroup))
    nodes = tree.nodes
    root = tree.root
    weight = {leaf: 0 if leaf in outgroup else 1 + len(members.get(leaf, []))
              for leaf in tree.leaves}
    below: dict[str, int] = {}
    for node in reversed(nodes):
        below[node] = weight.get(node, 0) + sum(below[child] for child in tree.children[node])

    removed: dict[str, str] = {}
    kept: set[str] = {root}
    for node in nodes:
        if node == root or not tree.children[node]:
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

    # The deepest kept node above each node, found from the root down.
    kept_above: dict[str, str] = {}
    for node in nodes[1:]:
        parent = tree.parents[node]
        kept_above[node] = parent if parent in kept else kept_above[parent]

    parent_of = {node: (None if node == root else kept_above[node])
                 for node in nodes if node in kept}
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
        parent_of = {node: (None if parent is None else joined_into.get(parent, parent))
                     for node, parent in parent_of.items() if node not in joined_into}
        genomes = {genome: joined_into.get(node, node) for genome, node in genomes.items()}

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
        ordered = sorted(tree.children[node], key=lambda child: (
            -count[child], child not in first, first.get(child, ""), child))
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
                    prefix: str) -> list[str]:
    """Write the hierarchy, the genome map and a node table."""
    os.makedirs(outdir, exist_ok=True)
    hierarchy = collapsed.hierarchy
    count = collapsed.genome_count()
    paths = [os.path.join(outdir, f"{prefix}_{name}")
             for name in ("hierarchy.tsv", "genomes.tsv", "nodes.tsv")]

    with open(paths[0], "w", newline="\n") as handle:
        handle.write("child\tparent\n")
        handle.writelines(f"{node}\t{hierarchy.parents[node] or ''}\n" for node in hierarchy.nodes)

    with open(paths[1], "w", newline="\n") as handle:
        handle.writelines(f"{genome}\t{node}\n" for genome, node in sorted(collapsed.genomes.items(), key=lambda item: (item[1], item[0])))

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
    clade: str | None = None
    jaccard: float = 0.0
    carriers: int = 0
    shared: int = 0
    clade_only: int = 0
    carriers_only: int = 0
    genomes: frozenset[str] = frozenset()
    reason: str = ""
    outside: bool = False
    """Whether the match is the genomes outside ``clade`` rather than inside it."""
    group: frozenset[str] = frozenset()
    """The genomes the node was matched to."""
    uncalled: frozenset[str] = frozenset()
    """Genomes of the tree without a base at the SNP the match rests on, which
    therefore say nothing about whether they belong to the node."""
    rank: int = 1
    """The rank of the node's canSNP the match rests on."""

    @property
    def described(self) -> str:
        return f"the genomes outside {self.clade}" if self.outside else str(self.clade)

    HEADER = ("node", "tree_node", "side", "jaccard", "carriers", "shared", "clade_only",
              "carriers_only", "genomes_placed", "note")


@dataclass
class Placement:
    """Where the genomes sit in an existing scheme, read from a genome tree."""

    clusters: dict[str, list[str]]
    """The genomes placed at each node of the scheme, by alignment record name."""
    matches: list[Match]
    findings: list[Finding]


def _clade(tree: Hierarchy, node: str, records: Mapping[str, str]) -> frozenset[str]:
    """The alignment records at or below a node of the tree."""
    return frozenset(records[leaf] for leaf in tree.leaves_below(node) if leaf in records)


def _clade_sizes(tree: Hierarchy, records: Mapping[str, str]) -> dict[str, int]:
    """The number of alignment records at or below each node of the tree.

    Counted rather than collected, since the sets of every node of a deep
    tree together grow with the square of its size.
    """
    sizes: dict[str, int] = {}
    for node in reversed(tree.nodes):
        sizes[node] = (node in records) + sum(sizes[child] for child in tree.children[node])
    return sizes


def place(scheme: Scheme, tree: Hierarchy, alignment: Mapping[str, str], reference: Reference,
          genome: str | None = None, min_match: float = 0.5,
          members: Mapping[str, Sequence[str]] | None = None) -> Placement:
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
    findings: list[Finding] = []
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

    member_record: dict[str, str] = {}
    folded = {leaf: names for leaf, names in (members or {}).items() if leaf in leaf_record}
    # Matched against the records once, not leaf by leaf.
    found, _missing, _ambiguous = resolve_names(
        list(dict.fromkeys(name for names in folded.values() for name in names)), records)
    for leaf, names in folded.items():
        for name in names:
            if name in found:
                member_record[found[name]] = leaf

    nodes = tree.nodes
    root = tree.root
    sizes = _clade_sizes(tree, leaf_record)
    in_tree = frozenset(leaf_record[leaf] for leaf in tree.leaves if leaf in leaf_record)

    matches: dict[str, Match] = {}
    annotated = sorted({snp.genome for snp in scheme.snps})
    if genome is None and len(annotated) > 1:
        raise ValueError(
            f"the scheme annotates {len(annotated)} reference genomes ({', '.join(annotated)}), "
            f"so --genome must say which one the alignment is against")
    if genome is not None and genome not in annotated:
        raise ValueError(
            f"the scheme annotates no SNP on genome '{genome}'; it names "
            f"{', '.join(annotated) or 'none'}")
    # Every canSNP of a node is tried, best rank first, and the one that
    # matches a side of a branch best is kept: one homoplastic or uncalled
    # canSNP then does not misplace a node whose other canSNPs mark it.
    snps: dict[str, list] = {}
    for snp in sorted(scheme.checked_snps, key=lambda snp: (snp.rank, snp.line)):
        if genome is not None and snp.genome != genome:
            continue
        snps.setdefault(snp.node, []).append(snp)

    def match_snp(node: str, snp) -> Match:
        if not 1 <= snp.position <= reference.length:
            return Match(node, reason="the position is outside the reference")
        column = mapping.to_column(snp.position - 1)
        if column is None:
            return Match(node, reason="the position is not aligned")

        called: set[str] = set()
        carriers: set[str] = set()
        for record in in_tree:
            base = alignment[record][column].upper()
            if base in BASES:
                called.add(record)
                if base == snp.derived.upper():
                    carriers.add(record)
        if not carriers:
            return Match(node, reason="no genome of the tree carries the derived base")

        # Count, for every node of the tree, the genomes below it that have a
        # base here and those that carry the derived one. Each branch splits
        # the genomes into the clade below it and the rest; both sides are
        # candidates, so the match does not depend on where the tree is rooted.
        n_called = len(called)
        n_carriers = len(carriers)
        below_called: dict[str, int] = {}
        below_carriers: dict[str, int] = {}
        for tree_node in reversed(nodes):
            record = leaf_record.get(tree_node)
            own_called = 1 if record in called else 0
            own_carrier = 1 if record in carriers else 0
            below_called[tree_node] = own_called + sum(
                below_called[child] for child in tree.children[tree_node])
            below_carriers[tree_node] = own_carrier + sum(
                below_carriers[child] for child in tree.children[tree_node])

        best: tuple[float, int, int, str] | None = None
        for tree_node in nodes:
            inside_called, inside_carriers = below_called[tree_node], below_carriers[tree_node]
            sides = [(inside_called, inside_carriers, 0)]
            # The genomes outside a clade are the other side of its branch.
            # Those outside one leaf are a side every tree has, which says
            # nothing about clades, and on an unresolved tree would stand in
            # for one; they are taken only as an exact match, as when the tree
            # is unrooted and one outgroup genome lies outside the scheme's
            # root, and rank after a clade that matches as well.
            if tree_node != root:
                sides.append((n_called - inside_called, n_carriers - inside_carriers,
                              1 if tree.children[tree_node] else 2))
            for side_called, shared, outside in sides:
                if not shared:
                    continue
                if outside == 2 and not (shared == n_carriers == side_called):
                    continue
                index = shared / (side_called + n_carriers - shared)
                size = sizes[tree_node] if not outside else len(in_tree) - sizes[tree_node]
                # The higher index first; then the side that is a clade of the
                # tree as rooted; then the smaller side.
                key = (-index, outside, size, tree_node)
                if best is None or key < best:
                    best = key

        index, outside, _size, clade = -best[0], bool(best[1]), best[2], best[3]
        inside = _clade(tree, clade, leaf_record)
        group = in_tree - inside if outside else inside
        genomes = group & called
        match = Match(node, clade, index, len(carriers), len(genomes & carriers),
                      len(genomes - carriers), len(carriers - genomes),
                      outside=outside, group=group, uncalled=frozenset(in_tree - called),
                      rank=snp.rank)
        if index < min_match:
            match.reason = f"the best clade matches with a Jaccard index of {index:.2f} only"
            match.clade = None
        elif outside:
            match.reason = (f"matches the genomes outside the tree's clade {clade}, so the tree "
                            f"is rooted inside this node's clade")
        return match

    for node in scheme.hierarchy.nodes:
        if node not in snps:
            matches[node] = Match(node, reason="no SNP for this node")
            continue
        tried = [match_snp(node, snp) for snp in snps[node]]
        # A matched clade before none, then the higher index, then the better rank.
        matches[node] = min(tried, key=lambda match: (match.clade is None, -match.jaccard,
                                                      match.rank))

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
        # Genomes without a base at the parent's SNP may belong to it: only a
        # genome called there and left outside contradicts the scheme.
        if child_match.clade and parent_match and parent_match.clade \
                and not child_match.group <= parent_match.group | parent_match.uncalled:
            findings.append(problem(
                MODIFY, "placement-order", node,
                f"the node matches {child_match.described}, which is not inside "
                f"{parent_match.described}, which its parent {parent} matches. The tree and the "
                f"scheme disagree on where this branch belongs"))

    depth = {node: len(scheme.hierarchy.path(node)) for node in scheme.hierarchy.nodes}
    matched_nodes = [node for node, match in matches.items() if match.clade]
    clusters: dict[str, list[str]] = {}
    split: list[str] = []
    placed_at: dict[str, str] = {}
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


def place_by_clade(tree: Hierarchy, scheme_genomes: Mapping[str, str],
                   parents: Mapping[str, str | None],
                   outgroups: Sequence[str]) -> tuple[dict[str, str], list[str]]:
    """Place the leaves of a tree at the deepest scheme node whose clade holds them.

    Where ``place`` reads a node's clade from where its canSNP's derived base
    lies, this reads it from the genomes the scheme was made from: the clade
    of a scheme node is that of the most recent common ancestor, in the tree
    rooted on the outgroups, of the scheme's genomes at or below the node. A
    genome added to the tree, such as one held out when the scheme was made,
    is then placed at the deepest scheme node whose clade holds it, without
    reference to any canSNP. The clade of a node with a single genome in the
    tree is that genome's leaf alone, so no other genome is placed there.

    Args:
        tree: A genome tree holding the scheme's genomes, the genomes to place
            and the outgroups.
        scheme_genomes: The node of each of the scheme's genomes, named by
            leaf or as ``resolve_names`` matches a leaf, such as by accession.
        parents: The parent of each node of the scheme, None for its root.
        outgroups: The leaves to root the tree on, as ``reroot`` takes them.

    Returns:
        The node of every leaf that is neither a scheme genome nor an
        outgroup and lies in a clade; and, sorted, those that lie in none.
        Among scheme nodes of the same depth whose clades both hold a leaf,
        which happens only where the tree does not keep a node's genomes
        together, the one with the smaller clade is taken.

    Raises:
        HierarchyError: If an outgroup is not a leaf, or the outgroups do not
            form a clade of the tree.
        ValueError: If a scheme genome sits at a node the scheme lacks.
    """
    unknown = sorted({node for node in scheme_genomes.values() if node not in parents})
    if unknown:
        raise ValueError(
            f"{len(unknown)} nodes of the scheme genomes are not nodes of the scheme: "
            f"{', '.join(unknown[:5])}. A scheme exported with --prune lacks the nodes it left "
            f"out; its genomes file places the genomes at the nodes it kept")
    outgroup = set(outgroup_leaves(tree, outgroups))
    tree = reroot(tree, sorted(outgroup))
    matched, unmatched, ambiguous = resolve_names(list(scheme_genomes), tree.leaves)
    if unmatched:
        logging.info(f"{len(unmatched)} scheme genomes are not in the tree: "
                     f"{', '.join(sorted(unmatched)[:5])}")
    if ambiguous:
        logging.warning(f"{len(ambiguous)} scheme genomes match more than one leaf of the tree "
                        f"and were left out: {', '.join(sorted(ambiguous)[:5])}")
    leaf_node = {leaf: scheme_genomes[genome] for genome, leaf in matched.items()
                 if leaf not in outgroup}

    # The scheme's genomes in the tree at or below each scheme node.
    below: dict[str, set[str]] = {}
    for leaf, node in leaf_node.items():
        walk: str | None = node
        while walk is not None:
            below.setdefault(walk, set()).add(leaf)
            walk = parents.get(walk)
    depth: dict[str, int] = {}
    for node in below:
        walk, steps = node, 0
        while walk is not None:
            steps += 1
            walk = parents.get(walk)
        depth[node] = steps

    # The most recent common ancestor of a set of leaves is the deepest node
    # of the tree that holds them all; nodes come parents first, so the last
    # such node in that order is the deepest.
    nodes = tree.nodes
    clade_of: dict[str, frozenset[str]] = {}
    for scheme_node, leaves in below.items():
        count: dict[str, int] = {}
        for node in reversed(nodes):
            count[node] = (node in leaves) + sum(count[child] for child in tree.children[node])
        mrca = next(node for node in reversed(nodes) if count[node] == len(leaves))
        clade_of[scheme_node] = frozenset(tree.leaves_below(mrca))

    order = sorted(clade_of, key=lambda node: (-depth[node], len(clade_of[node]), node))
    placed: dict[str, str] = {}
    unplaced: list[str] = []
    for leaf in tree.leaves:
        if leaf in leaf_node or leaf in outgroup:
            continue
        node = next((node for node in order if leaf in clade_of[node]), None)
        if node is None:
            unplaced.append(leaf)
        else:
            placed[leaf] = node
    logging.info(
        f"Found the clades of {len(clade_of)} scheme nodes from {len(leaf_node)} scheme genomes "
        f"in the tree, and placed {len(placed)} of {len(placed) + len(unplaced)} other genomes")
    return placed, sorted(unplaced)


def write_placement(placement: Placement, genome_map: str, report: str | None = None) -> None:
    """Write the genome map, and a report of how each node was matched."""
    with open(genome_map, "w", newline="\n") as handle:
        for node, genomes in sorted(placement.clusters.items()):
            handle.writelines(f"{genome}\t{node}\n" for genome in sorted(genomes))
    if report:
        with open(report, "w", newline="\n") as handle:
            handle.write("\t".join(Match.HEADER) + "\n")
            handle.writelines("\t".join([
                    match.node, match.clade or "",
                    "" if not match.clade else ("outside" if match.outside else "inside"),
                    f"{match.jaccard:.3f}", str(match.carriers),
                    str(match.shared), str(match.clade_only), str(match.carriers_only),
                    str(len(match.genomes)), match.reason]) + "\n" for match in placement.matches)

