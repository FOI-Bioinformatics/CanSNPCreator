"""Placing genomes on an existing scheme.

To check an existing scheme against genomes, each genome must sit at a node of
that scheme. Typing the genomes with the scheme itself places them
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

import logging
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from cansnpscorer.findings import CHECK, MODIFY, Finding, problem, recommendation
from cansnpscorer.helpers.names import resolve_names
from cansnpscorer.hierarchy import Hierarchy
from cansnpscorer.model import AlignmentMap, Reference
from cansnpscorer.scheme import Scheme
from cansnpscorer.trees.rooting import outgroup_leaves, reroot

BASES = frozenset("ACGT")


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

    HEADER = (
        "node",
        "tree_node",
        "side",
        "jaccard",
        "carriers",
        "shared",
        "clade_only",
        "carriers_only",
        "genomes_placed",
        "note",
    )


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


def place(
    scheme: Scheme,
    tree: Hierarchy,
    alignment: Mapping[str, str],
    reference: Reference,
    genome: str | None = None,
    min_match: float = 0.5,
    members: Mapping[str, Sequence[str]] | None = None,
) -> Placement:
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
        return Placement(
            {},
            [],
            [problem(CHECK, "placement", "", "the tree file of the scheme is not a single rooted tree")],
        )

    mapping = AlignmentMap(next(iter(alignment.values())), reference)
    records = list(alignment)
    matched, unmatched, ambiguous = resolve_names(tree.leaves, records)
    if unmatched or ambiguous:
        findings.append(
            recommendation(
                CHECK,
                "tree-genomes",
                "",
                f"{len(unmatched) + len(ambiguous)} leaves of the tree match no single genome of "
                f"the alignment and are left out: {', '.join(sorted(unmatched + ambiguous)[:5])}",
            )
        )
    leaf_record = dict(matched)

    member_record: dict[str, str] = {}
    folded = {leaf: names for leaf, names in (members or {}).items() if leaf in leaf_record}
    # Matched against the records once, not leaf by leaf.
    found, _missing, _ambiguous = resolve_names(
        list(dict.fromkeys(name for names in folded.values() for name in names)), records
    )
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
            f"so --genome must say which one the alignment is against"
        )
    if genome is not None and genome not in annotated:
        raise ValueError(
            f"the scheme annotates no SNP on genome '{genome}'; it names {', '.join(annotated) or 'none'}"
        )
    # Every canSNP of a node is tried, best rank first, and the one that
    # matches a side of a branch best is kept: one homoplastic or uncalled
    # canSNP then does not misplace a node whose other canSNPs mark it.
    snps: dict[str, list] = {}
    for snp in sorted(scheme.checked_snps, key=lambda snp: (snp.rank, snp.line)):
        if genome is not None and snp.genome != genome:
            continue
        snps.setdefault(snp.node, []).append(snp)

    def match_snp(node: str, snp) -> Match:
        if snp.sequence not in reference.names:
            return Match(node, reason=f"the reference holds no sequence '{snp.sequence}'")
        if not 1 <= snp.position <= reference.contig(snp.sequence).length:
            return Match(node, reason="the position is outside its reference sequence")
        column = mapping.column_at(snp.sequence, snp.position)
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
            node_record = leaf_record.get(tree_node)
            own_called = 1 if node_record in called else 0
            own_carrier = 1 if node_record in carriers else 0
            below_called[tree_node] = own_called + sum(
                below_called[child] for child in tree.children[tree_node]
            )
            below_carriers[tree_node] = own_carrier + sum(
                below_carriers[child] for child in tree.children[tree_node]
            )

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
                sides.append(
                    (
                        n_called - inside_called,
                        n_carriers - inside_carriers,
                        1 if tree.children[tree_node] else 2,
                    )
                )
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

        # The root, taken as the clade below it, holds every carrier, and
        # there is at least one, so a side was kept.
        assert best is not None
        index, outside, _size, clade = -best[0], bool(best[1]), best[2], best[3]
        inside = _clade(tree, clade, leaf_record)
        group = in_tree - inside if outside else inside
        genomes = group & called
        match = Match(
            node,
            clade,
            index,
            len(carriers),
            len(genomes & carriers),
            len(genomes - carriers),
            len(carriers - genomes),
            outside=outside,
            group=group,
            uncalled=frozenset(in_tree - called),
            rank=snp.rank,
        )
        if index < min_match:
            match.reason = f"the best clade matches with a Jaccard index of {index:.2f} only"
            match.clade = None
        elif outside:
            match.reason = (
                f"matches the genomes outside the tree's clade {clade}, so the tree "
                f"is rooted inside this node's clade"
            )
        return match

    for node in scheme.hierarchy.nodes:
        if node not in snps:
            matches[node] = Match(node, reason="no SNP for this node")
            continue
        tried = [match_snp(node, snp) for snp in snps[node]]
        # A matched clade before none, then the higher index, then the better rank.
        matches[node] = min(tried, key=lambda match: (match.clade is None, -match.jaccard, match.rank))

    for node, match in matches.items():
        if match.clade is None:
            findings.append(
                recommendation(
                    CHECK,
                    "unplaced-node",
                    node,
                    f"the node could not be matched to a clade of the tree: {match.reason}. No "
                    f"genome is placed at it",
                )
            )
        elif match.jaccard < 1.0:
            findings.append(
                recommendation(
                    CHECK,
                    "partial-match",
                    node,
                    f"the node matches {match.described} with a Jaccard index of "
                    f"{match.jaccard:.2f}: {match.clade_only} genomes of the clade lack the derived "
                    f"base, and {match.carriers_only} genomes outside it carry it",
                )
            )

    rerooted = sorted(node for node, match in matches.items() if match.clade and match.outside)
    if rerooted:
        findings.append(
            recommendation(
                CHECK,
                "tree-root",
                ", ".join(rerooted[:5]),
                f"{len(rerooted)} nodes match the genomes outside a clade of the tree rather than a "
                f"clade, so the tree is rooted inside them. The placement does not depend on it, "
                f"but a tree rooted on an outgroup avoids the question",
            )
        )

    for node in scheme.hierarchy.nodes:
        parent = scheme.hierarchy.parents[node]
        child_match, parent_match = matches[node], matches.get(parent) if parent else None
        # Genomes without a base at the parent's SNP may belong to it: only a
        # genome called there and left outside contradicts the scheme.
        if (
            child_match.clade
            and parent_match
            and parent_match.clade
            and not child_match.group <= parent_match.group | parent_match.uncalled
        ):
            findings.append(
                problem(
                    MODIFY,
                    "placement-order",
                    node,
                    f"the node matches {child_match.described}, which is not inside "
                    f"{parent_match.described}, which its parent {parent} matches. The tree and the "
                    f"scheme disagree on where this branch belongs",
                )
            )

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
        findings.append(
            problem(
                CHECK,
                "placement-conflict",
                "",
                f"{len(split)} genomes lie in the clades of nodes on different branches of the "
                f"scheme, and were placed at the deepest of them: {', '.join(split[:5])}",
            )
        )

    for match in matches.values():
        match.genomes = frozenset(clusters.get(match.node, []))

    logging.info(
        f"Matched {len(matched_nodes)} of {len(matches)} nodes of the scheme to clades of the "
        f"tree, and placed {len(placed_at)} of {len(in_tree)} genomes"
    )
    return Placement(clusters, list(matches.values()), findings)


def place_by_clade(
    tree: Hierarchy,
    scheme_genomes: Mapping[str, str],
    parents: Mapping[str, str | None],
    outgroups: Sequence[str],
) -> tuple[dict[str, str], list[str]]:
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
            f"out; its genomes file places the genomes at the nodes it kept"
        )
    outgroup = set(outgroup_leaves(tree, outgroups))
    tree = reroot(tree, sorted(outgroup))
    matched, unmatched, ambiguous = resolve_names(list(scheme_genomes), tree.leaves)
    if unmatched:
        logging.info(
            f"{len(unmatched)} scheme genomes are not in the tree: {', '.join(sorted(unmatched)[:5])}"
        )
    if ambiguous:
        logging.warning(
            f"{len(ambiguous)} scheme genomes match more than one leaf of the tree "
            f"and were left out: {', '.join(sorted(ambiguous)[:5])}"
        )
    leaf_node = {leaf: scheme_genomes[genome] for genome, leaf in matched.items() if leaf not in outgroup}

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
        owner = next((node for node in order if leaf in clade_of[node]), None)
        if owner is None:
            unplaced.append(leaf)
        else:
            placed[leaf] = owner
    logging.info(
        f"Found the clades of {len(clade_of)} scheme nodes from {len(leaf_node)} scheme genomes "
        f"in the tree, and placed {len(placed)} of {len(placed) + len(unplaced)} other genomes"
    )
    return placed, sorted(unplaced)


def write_placement(placement: Placement, genome_map: str, report: str | None = None) -> None:
    """Write the genome map, and a report of how each node was matched."""
    with open(genome_map, "w", newline="\n") as handle:
        for node, genomes in sorted(placement.clusters.items()):
            handle.writelines(f"{genome}\t{node}\n" for genome in sorted(genomes))
    if report:
        with open(report, "w", newline="\n") as handle:
            handle.write("\t".join(Match.HEADER) + "\n")
            handle.writelines(
                "\t".join(
                    [
                        match.node,
                        match.clade or "",
                        "" if not match.clade else ("outside" if match.outside else "inside"),
                        f"{match.jaccard:.3f}",
                        str(match.carriers),
                        str(match.shared),
                        str(match.clade_only),
                        str(match.carriers_only),
                        str(len(match.genomes)),
                        match.reason,
                    ]
                )
                + "\n"
                for match in placement.matches
            )
