"""The part of an audit that needs the genomes.

Checking a scheme against its reference genome finds bookkeeping errors. It
cannot find the property that defines a canSNP: that the substitution is
carried by every isolate below its node and by no other genome. That needs the
base of every genome at every SNP, which a whole genome alignment holds.

Given the alignment and a genome map placing its genomes at nodes of the
scheme, this module answers three questions:

1. Is each SNP of the scheme canonical in these genomes? It reports isolates
   below the node that do not carry the derived base, genomes outside the node
   that do, SNPs whose two alleles appear swapped, and genomes with no base at
   the site.
2. For a node whose SNP fails, is there a canonical one to replace it with?
   And is there any site at all that separates exactly the node's isolates, or
   do these genomes not form that clade?
3. Should the scheme be extended? A group of genomes that shares substitutions
   no other genome carries, and that the scheme has no node for, is a clade the
   scheme does not resolve. Sites that fit no grouping compatible with the
   scheme are counted too, since a genome placed at the wrong node, or a
   recombinant region, shows itself as many such sites.

Genomes the genome map places nowhere are counted as outside every node, as
the extract subcommand counts them, and are reported once.
"""

import logging
from collections import Counter
from collections.abc import Mapping, Sequence

from cansnpscorer.findings import CHECK, EXTEND, MODIFY, Finding, problem, recommendation
from cansnpscorer.helpers.data_utils import snp_extractor
from cansnpscorer.helpers.names import resolve_clusters
from cansnpscorer.model import AlignmentMap, Reference, ReferenceError
from cansnpscorer.scheme import Scheme
from cansnpscorer.typing_report import ANCESTRAL, DERIVED, derived_off, set_path

BASES = frozenset("ACGT")

#: At most this many genome names are written into one finding.
SHOWN = 5
#: A canSNP without a base in more than this share of the genomes is reported
#: on its own; the others are counted in one finding.
MISSING_SHARE = 0.05
#: At most this many new clades are reported one by one; the rest are counted.
REPORTED_CLADES = 50


def _names(genomes: Sequence[str]) -> str:
    shown = sorted(genomes)[:SHOWN]
    more = len(genomes) - len(shown)
    return ", ".join(shown) + (f" and {more} more" if more > 0 else "")


def _where(reference: Reference, joined_position: int) -> str:
    """A 0-based position along the joined reference sequences, as sequence:position."""
    contig, local = reference.locate(joined_position)
    return f"{contig.name}:{local + 1}"


class Placement:
    """Where the genome map puts each genome, and the isolates below each node."""

    def __init__(self, scheme: Scheme, clusters: Mapping[str, Sequence[str]],
                 genomes: Sequence[str]) -> None:
        hierarchy = scheme.hierarchy
        self.hierarchy = hierarchy
        self.unknown_clusters = sorted(set(clusters) - set(hierarchy.parents))
        # A genome the map puts at a node the scheme does not hold, such as
        # one export --prune left out, has no known place in this scheme. It
        # is neither below any node nor known to be outside one, so it tests
        # no SNP: counted as outside, it would carry the derived base of
        # every node above the one it was placed at, and each of those SNPs
        # would be reported as homoplastic.
        self.unplaced: list[str] = sorted(
            {name for node in self.unknown_clusters for name in clusters[node]} & set(genomes))
        without_place = set(self.unplaced)
        self.genomes: list[str] = [name for name in genomes if name not in without_place]
        known = set(self.genomes)

        self.missing_genomes = sorted(
            {name for node, names in clusters.items() for name in names} - set(genomes))

        usable = {
            node: [name for name in names if name in known]
            for node, names in clusters.items() if node in hierarchy
        }
        self.node_of: dict[str, str] = {}
        for node, names in usable.items():
            for name in names:
                self.node_of.setdefault(name, node)

        below = hierarchy.isolates(usable)
        self.below: dict[str, frozenset[str]] = {
            node: frozenset(names) for node, names in below.items()}
        self.ungrouped = sorted(known - set(self.node_of))
        self.node_by_set: dict[frozenset[str], str] = {}
        # Parents come before children, so the deepest node of a set is kept;
        # a node with one child holds the same isolates as that child.
        for node in hierarchy.nodes:
            if self.below[node]:
                self.node_by_set[self.below[node]] = node

    def lowest_holding(self, carriers: frozenset[str]) -> str | None:
        """The deepest node whose isolates include all of these genomes."""
        first = next(iter(carriers))
        node: str | None = self.node_of.get(first)
        while node is not None and not carriers <= self.below[node]:
            node = self.hierarchy.parents[node]
        return node

    def fits(self, carriers: frozenset[str], node: str) -> bool:
        """Whether a group of genomes below ``node`` respects the scheme's children.

        Each child of the node must lie wholly inside the group or wholly
        outside it. A group that cuts a child in two is not a clade of this
        scheme.
        """
        for child in self.hierarchy.children[node]:
            inside = self.below[child] & carriers
            if inside and inside != self.below[child]:
                return False
        return True


def _bases_at(alignment: Mapping[str, str], column: int,
              genomes: Sequence[str]) -> dict[str, list[str]]:
    by_base: dict[str, list[str]] = {}
    for genome in genomes:
        base = alignment[genome][column].upper()
        by_base.setdefault(base if base in BASES else "missing", []).append(genome)
    return by_base


def check_scheme_snps(scheme: Scheme, alignment: Mapping[str, str], mapping: AlignmentMap,
                      placement: Placement, genome: str | None) -> tuple[list[Finding],
                                                                           set[str]]:
    """Whether each SNP of the scheme is canonical in these genomes.

    Returns the findings and the nodes whose SNP is not canonical.
    """
    findings: list[Finding] = []
    failed: set[str] = set()
    untested: list[str] = []
    # Genomes without a base at a canSNP, summed over the canSNPs that lack
    # few of them: one finding per canSNP hid the rest of the report, and
    # which genomes are behind them was not visible at all.
    rarely_missing: list[tuple[str, int]] = []
    missing_by_genome: Counter = Counter()
    canonical = 0
    not_canonical = 0

    for snp in scheme.checked_snps:
        if genome is not None and snp.genome != genome:
            continue
        if snp.node not in placement.below:
            continue
        if (snp.sequence not in mapping.reference.names
                or not 1 <= snp.position <= mapping.reference.contig(snp.sequence).length):
            continue  # reported by the reference check
        inside = placement.below[snp.node]
        if not inside:
            if snp.typed:
                untested.append(snp.node)
            continue

        column = mapping.column_at(snp.sequence, snp.position)
        if column is None:
            findings.append(recommendation(
                CHECK, "not-aligned", snp.node, snp.tag(
                f"no column of the alignment covers {snp.sequence}:{snp.position}, so the SNP "
                f"cannot be checked in these genomes")))
            continue

        outside = [name for name in placement.genomes if name not in inside]
        if not outside:
            # Every genome sits below this node, which is the root in practice.
            # A SNP all of them carry separates nothing in this alignment, and
            # its polarity cannot be told without a genome outside the node.
            findings.append(recommendation(
                CHECK, "no-outgroup", snp.node, snp.tag(
                f"every genome of the alignment sits below this node, so nothing tests whether "
                f"its SNP at {snp.sequence}:{snp.position} marks it, nor which base is "
                f"ancestral. Add an "
                f"outgroup genome to the alignment, placed at no node")))
            continue

        within = _bases_at(alignment, column, sorted(inside))
        beyond = _bases_at(alignment, column, outside)
        derived, ancestral = snp.derived.upper(), snp.ancestral.upper()
        where = f"at {snp.sequence}:{snp.position}"

        called_inside = {base for base in within if base != "missing"}
        bad = False
        not_derived = [name for base, names in within.items()
                       if base not in (derived, "missing") for name in names]
        derived_outside = beyond.get(derived, [])
        missing = within.get("missing", []) + beyond.get("missing", [])

        if called_inside == {ancestral} and derived_outside:
            bad = True
            findings.append(problem(
                MODIFY, "polarity", snp.node, snp.tag(
                f"the isolates below the node carry '{ancestral}' {where}, which the scheme gives "
                f"as the ancestral base, and {len(derived_outside)} genomes outside it carry the "
                f"derived '{derived}'. The two columns appear swapped")))
            not_canonical += 1
            if snp.typed:
                failed.add(snp.node)
            continue

        if not_derived:
            bad = True
            findings.append(problem(
                MODIFY, "not-canonical", snp.node, snp.tag(
                f"{len(not_derived)} of {len(inside)} isolates below the node do not carry the "
                f"derived '{derived}' {where}: {_names(not_derived)}. CanSNPer2 would place "
                f"them above this node")))

        if derived_outside:
            bad = True
            placed = Counter(placement.node_of.get(name, "no node") for name in derived_outside)
            nodes = ", ".join(f"{node} ({count})" for node, count in placed.most_common(SHOWN))
            findings.append(problem(
                MODIFY, "homoplasy", snp.node, snp.tag(
                f"{len(derived_outside)} genomes outside the node carry the derived "
                f"'{derived}' {where}: {_names(derived_outside)}; they are placed at {nodes}. "
                f"The SNP does not mark this branch alone")))

        for name in missing:
            missing_by_genome[name] += 1
        if missing and len(missing) <= MISSING_SHARE * len(placement.genomes):
            rarely_missing.append((snp.node, len(missing)))
        elif missing:
            findings.append(recommendation(
                CHECK, "missing-data", snp.node, snp.tag(
                f"{len(missing)} genomes have no base {where}: {_names(missing)}. The SNP cannot "
                f"be called in them, which is how an accessory or poorly aligned region shows "
                f"itself")))

        if not bad:
            canonical += 1
            continue
        not_canonical += 1
        if snp.typed:
            # Only the canSNP CanSNPer2 calls is offered a replacement.
            failed.add(snp.node)

    if rarely_missing:
        rarely_missing.sort(key=lambda item: (-item[1], item[0]))
        nodes = ", ".join(f"{node} ({count})" for node, count in rarely_missing[:SHOWN])
        findings.append(recommendation(
            CHECK, "missing-data", "",
            f"{len(rarely_missing)} further canSNPs have no base in a few genomes each, at most "
            f"{MISSING_SHARE:.0%} of them; the most affected are at {nodes}"))
    if missing_by_genome:
        worst = ", ".join(f"{name} ({count})"
                          for name, count in sorted(missing_by_genome.items(),
                                                    key=lambda item: (-item[1], item[0]))[:SHOWN])
        findings.append(recommendation(
            CHECK, "missing-data-genomes", "",
            f"{len(missing_by_genome)} genomes have no base at one canSNP or more, and the ones "
            f"lacking most are {worst}. A genome that lacks many is more likely incomplete or "
            f"poorly aligned than evidence against the canSNPs"))

    if untested:
        findings.append(recommendation(
            CHECK, "untested", ", ".join(sorted(untested)[:SHOWN]),
            f"{len(untested)} nodes have no genome placed at or below them, so their SNPs were "
            f"not checked"))

    logging.info(
        f"{canonical} SNPs of the scheme are canonical in these genomes, and {not_canonical} "
        f"are not")
    return findings, failed


#: A genome that contradicts the scheme at this many nodes or more is reported.
GENOME_CONFLICTS = 3


def _lineage(hierarchy, node: str) -> list[str]:
    """The node and every node above it."""
    nodes: list[str] = []
    while node is not None:
        nodes.append(node)
        node = hierarchy.parents[node]
    return nodes


def check_genome_paths(scheme: Scheme, alignment: Mapping[str, str], mapping: AlignmentMap,
                       placement: Placement, genome: str | None,
                       min_conflicts: int = GENOME_CONFLICTS) -> list[Finding]:
    """Whether each genome's own canSNPs agree with each other and with its placement.

    The canSNPs of a node give its state in a genome: derived when more of
    them carry the derived base than the ancestral one, ancestral when fewer,
    unknown otherwise. The genome's path leads to the derived node with no
    node above ancestral that has the most derived nodes on its way from the
    root, as the typing report's set rule takes it; two such nodes on
    different branches leave the path at the nodes they share. A node of
    unknown state may lie on the path, as when its canSNPs have no base in
    the genome; CanSNPer2 2.3.0 passes such a node but calls neither it nor
    the node directly below it, so the path is where the genome's canSNPs
    place it, which can lie below where CanSNPer2 would type it.

    Two things are reported, as they explain failures check_scheme_snps
    attributes to SNPs:

    - genomes whose states are consistent but lead to another node than the
      genome map gives. A genome map that puts a whole RepGenR2 cluster at its
      representative's node does this when the cluster holds genomes from both
      sides of a node of the scheme. With one canSNP at the node they are
      typed at, the same pattern arises when that SNP is homoplastic, so the
      number of canSNPs the typing rests on is given;
    - genomes that contradict the scheme at ``min_conflicts`` nodes or more:
      derived at a node off their path, or ancestral on it. On 1,229 Brucella
      genomes, half of such contradictions came from 5% of the genomes, which
      points to base errors or mixed assemblies rather than to the SNPs.
    """
    hierarchy = placement.hierarchy
    by_node: dict[str, list[tuple[int, str, str]]] = {}
    for snp in scheme.checked_snps:
        if genome is not None and snp.genome != genome:
            continue
        if snp.node not in hierarchy.parents:
            continue
        column = mapping.column_at(snp.sequence, snp.position)
        if column is not None:
            by_node.setdefault(snp.node, []).append(
                (column, snp.derived.upper(), snp.ancestral.upper()))
    if not by_node:
        return []

    lineage = {node: _lineage(hierarchy, node) for node in hierarchy.parents}
    misplaced: dict[tuple[str, str], list[tuple[str, int]]] = {}
    discordant: list[tuple[str, int]] = []
    for name in placement.genomes:
        placed = placement.node_of.get(name)
        if placed is None:
            continue
        sequence = alignment[name]
        derived_at: dict[str, int] = {}
        ancestral_at: set[str] = set()
        for node, snps in by_node.items():
            derived = sum(1 for column, d, _a in snps if sequence[column].upper() == d)
            ancestral = sum(1 for column, _d, a in snps if sequence[column].upper() == a)
            if derived > ancestral:
                derived_at[node] = derived
            elif ancestral > derived:
                ancestral_at.add(node)

        states = {**dict.fromkeys(ancestral_at, ANCESTRAL), **dict.fromkeys(derived_at, DERIVED)}
        typed, path = set_path(states, lineage)
        conflicts = (derived_off(states, lineage, typed, path)
                     + sum(1 for node in ancestral_at if node in path))
        if conflicts >= min_conflicts:
            discordant.append((name, conflicts))
            continue
        if conflicts or typed is None:
            continue
        expected = set(lineage[placed])
        against = (any(node not in expected for node in derived_at)
                   or any(node in expected for node in ancestral_at))
        if against:
            misplaced.setdefault((placed, typed), []).append((name, derived_at[typed]))

    findings: list[Finding] = []
    for (placed, typed), genomes in sorted(misplaced.items()):
        names = [name for name, _count in genomes]
        least = min(count for _name, count in genomes)
        findings.append(recommendation(
            CHECK, "misplaced", placed,
            f"{len(genomes)} genomes placed at {placed} carry the canSNPs of {typed} and none that "
            f"contradicts it, so they would be typed at {typed}: {_names(names)}. Either the "
            f"genome map places them wrongly, as it does when a cluster of genomes spans a node "
            f"of the scheme, or the canSNPs that lead to {typed} are not canonical; the typing "
            f"rests on at least {least} canSNP{'s' if least != 1 else ''} of {typed}"))
    if discordant:
        discordant.sort(key=lambda item: (-item[1], item[0]))
        shown = ", ".join(f"{name} ({count})" for name, count in discordant[:SHOWN])
        more = len(discordant) - min(len(discordant), SHOWN)
        findings.append(recommendation(
            CHECK, "discordant-genome", "",
            f"{len(discordant)} genomes contradict the scheme at {min_conflicts} or more nodes, "
            f"carrying the derived base off their path or the ancestral one on it: {shown}"
            + (f" and {more} more" if more > 0 else "")
            + ". Base errors, contamination or a mixed assembly show themselves this way; such "
              "genomes are weak evidence against a SNP"))
    logging.info(
        f"{sum(len(found) for found in misplaced.values())} genomes would be typed at another "
        f"node than the genome map gives, and {len(discordant)} contradict the scheme at "
        f"{min_conflicts} or more nodes")
    return findings


def scan_alignment(alignment: Mapping[str, str], mapping: AlignmentMap,
                   placement: Placement) -> tuple[dict[str, list[tuple[int, str, str]]],
                                                  dict[tuple[str, frozenset[str]],
                                                       list[tuple[int, str, str]]],
                                                  Counter, Counter, int]:
    """Sort every variable site by what it says about the scheme.

    Only sites where every genome has one of two bases are used, as the
    extract subcommand uses them.

    Returns:
        - canonical sites per node, as (position, derived, ancestral); the
          position is 0-based along the joined reference sequences, for use
          within this module only, and is written out as sequence:position
        - sites per candidate new clade, keyed by (parent node, genomes)
        - conflicting sites per node
        - conflicting sites per genome
        - the number of sites read
    """
    canonical: dict[str, list[tuple[int, str, str]]] = {}
    clades: dict[tuple[str, frozenset[str]], list[tuple[int, str, str]]] = {}
    conflicts_by_node: Counter = Counter()
    conflicts_by_genome: Counter = Counter()
    grouped = frozenset(placement.node_of)

    ordered = {name: alignment[name] for name in placement.genomes}
    sites = snp_extractor(ordered)

    for site in sites:
        column = site["index"]
        position = mapping.to_reference(column)
        if position is None:
            continue
        carriers: dict[str, set[str]] = {}
        for name in placement.genomes:
            carriers.setdefault(alignment[name][column].upper(), set()).add(name)
        if len(carriers) != 2:
            continue
        (first, first_set), (second, second_set) = sorted(
            carriers.items(), key=lambda item: (len(item[1]), item[0]))

        sides = ((first, frozenset(first_set), second), (second, frozenset(second_set), first))

        # A side that is a node of the scheme settles the site, whichever side
        # it is. Without an outgroup the two sides of a split are equally
        # rooted, so the side holding fewer genomes is not necessarily the
        # derived one: a site that separates C.7 from its sister nodes under
        # the root is C.7's canSNP, although the sister nodes hold fewer.
        placed = False
        for derived, group, ancestral in sides:
            node = placement.node_by_set.get(group)
            if node is not None:
                canonical.setdefault(node, []).append((position, derived, ancestral))
                placed = True
                break

        for derived, group, ancestral in sides if not placed else ():
            if len(group) < 2 or not group <= grouped:
                continue
            parent = placement.lowest_holding(group)
            if parent is not None and placement.fits(group, parent):
                clades.setdefault((parent, group), []).append((position, derived, ancestral))
                placed = True
                break

        if not placed:
            minor = frozenset(first_set)
            if minor <= grouped and len(minor) >= 2:
                node = placement.lowest_holding(minor)
                conflicts_by_node[node or "the root"] += 1
                for name in minor:
                    conflicts_by_genome[name] += 1

    return canonical, clades, conflicts_by_node, conflicts_by_genome, len(sites)


def check_extension(scheme: Scheme, mapping: AlignmentMap, placement: Placement,
                    canonical: dict[str, list[tuple[int, str, str]]],
                    clades: dict[tuple[str, frozenset[str]], list[tuple[int, str, str]]],
                    conflicts_by_node: Counter, conflicts_by_genome: Counter, failed: set[str],
                    sites: int, min_clade_size: int = 2, min_clade_snps: int = 1
                    ) -> list[Finding]:
    """Replacements for failed SNPs, nodes these genomes do not support, and new clades."""
    findings: list[Finding] = []
    reference = mapping.reference

    def described(entries: Sequence[tuple[int, str, str]]) -> str:
        shown = [f"{_where(reference, position)} ({ancestral}>{derived})"
                 for position, derived, ancestral in entries[:3]]
        return "; ".join(shown) + (" and others" if len(entries) > 3 else "")

    for node in scheme.hierarchy.nodes:
        if not placement.below.get(node):
            continue
        # Sites are credited to the deepest of nodes holding the same genomes,
        # and separate this node's genomes as well.
        found = canonical.get(placement.node_by_set.get(placement.below[node], node), [])
        if node in failed and found:
            findings.append(recommendation(
                MODIFY, "replacement", node,
                f"{len(found)} sites separate exactly the isolates below this node and could "
                f"replace its SNP: {described(found)}"))
        elif node in failed:
            findings.append(problem(
                CHECK, "unsupported-node", node,
                "no site of the alignment separates exactly the isolates below this node, so in "
                "these genomes it is not a clade. Check the placement of its isolates, or "
                "whether the node should be merged with its parent"))
        elif len(found) == 1 and node != scheme.hierarchy.root:
            # found[0] is 0-based along the joined reference sequences; the
            # scheme gives a 1-based position on one sequence.
            contig, local = reference.locate(found[0][0])
            site = (contig.name, local + 1)
            own = {(snp.sequence, snp.position) for snp in scheme.checked_snps
                   if snp.node == node}
            others = len(own - {site})
            if site in own:
                which = "one of the scheme's own canSNPs"
            else:
                which = f"which is not a canSNP of the scheme ({described(found)})"
            findings.append(recommendation(
                CHECK, "single-site", node,
                f"one site of the alignment at which every genome has a base separates exactly "
                f"the isolates below this node, {which}. "
                + (f"The node's {others} other canSNP{'s' if others != 1 else ''} lack a base "
                   f"in some genome or do not separate these isolates exactly, so the branch "
                   f"rests on little"
                   if others else
                   "A sequencing error in that one base moves a sample off the branch, and "
                   "nothing else in these genomes supports it")))

    ranked = sorted(clades.items(), key=lambda item: (-len(item[1]), item[0][0],
                                                      sorted(item[0][1])))

    # Each grouping fits the scheme on its own, but two of them can still
    # contradict each other: genomes shared without one holding the other. No
    # tree holds both, so the better supported one is kept and the other is
    # counted with the conflicting sites, as homoplasy or noise.
    accepted: list[frozenset[str]] = []
    candidates = []
    contradicted = 0
    for (parent, group), entries in ranked:
        if any(group & kept and not (group <= kept or kept <= group) for kept in accepted):
            contradicted += len(entries)
            conflicts_by_node[parent] += len(entries)
            for name in group:
                conflicts_by_genome[name] += len(entries)
            continue
        accepted.append(group)
        if len(group) >= min_clade_size and len(entries) >= min_clade_snps:
            candidates.append((parent, group, entries))

    for parent, group, entries in candidates[:REPORTED_CLADES]:
        children = [child for child in scheme.hierarchy.children[parent]
                    if placement.below[child] and placement.below[child] <= group]
        if children:
            shape = f"they span the child nodes {', '.join(children)} without the others"
        else:
            shape = f"they are {len(group)} of the {len(placement.below[parent])} isolates there"
        findings.append(recommendation(
            EXTEND, "new-clade", parent,
            f"{len(group)} genomes below this node share {len(entries)} SNPs no other genome "
            f"carries, and the scheme has no node for them; {shape}: {_names(sorted(group))}. "
            f"First sites: {described(entries)}"))

    if len(candidates) > REPORTED_CLADES:
        findings.append(recommendation(
            EXTEND, "new-clade", "",
            f"{len(candidates) - REPORTED_CLADES} further groups of genomes share SNPs that no "
            f"node of the scheme accounts for; raise --min-clade-snps to see the best "
            f"supported only"))

    conflicting = sum(conflicts_by_node.values())
    if conflicting:
        nodes = ", ".join(f"{node} ({count})" for node, count in conflicts_by_node.most_common(3))
        genomes = ", ".join(
            f"{name} ({count})" for name, count in conflicts_by_genome.most_common(3))
        contradicting = (
            f", {contradicted} of them because they group genomes against a better supported "
            f"grouping" if contradicted else "")
        findings.append(recommendation(
            CHECK, "conflicting-sites", "",
            f"{conflicting} of {sites} sites group the genomes in a way the scheme cannot "
            f"hold{contradicting}. Most fall below {nodes}, and the genomes most often involved "
            f"are {genomes}. A genome that dominates this list is likely placed at the wrong "
            f"node; many sites close together point to recombination"))

    return findings


def check_against_alignment(scheme: Scheme, alignment: Mapping[str, str], reference: Reference,
                            clusters: Mapping[str, Sequence[str]],
                            genome: str | None = None, min_clade_size: int = 2,
                            min_clade_snps: int = 1,
                            genome_conflicts: int = GENOME_CONFLICTS,
                            new_clades: bool = True) -> list[Finding]:
    """Every check that needs the genomes, as one list of findings.

    Args:
        scheme: The scheme as read from its source files.
        alignment: The whole genome alignment, reference first, in either of
            the layouts ``AlignmentMap`` accepts.
        reference: The reference genome the scheme is annotated against.
        clusters: The genomes placed at each node, as a genome map gives them.
        genome: Which of the scheme's reference genomes the alignment is
            against, when it names more than one.
        min_clade_size: The fewest genomes a new clade is reported for.
        min_clade_snps: The fewest SNPs a new clade is reported with.
        genome_conflicts: The fewest nodes a genome contradicts the scheme at
            for it to be reported.
        new_clades: Whether to look for clades the scheme does not resolve.
            That reads every column of the alignment; without it only the
            columns of the scheme's SNPs are read.
    """
    findings: list[Finding] = []

    if scheme.hierarchy is None:
        return [problem(CHECK, "alignment", "",
                        "the tree file is not a single rooted tree, so the genomes cannot be "
                        "compared against it")]

    if genome is None:
        genomes_annotated = sorted({snp.genome for snp in scheme.snps})
        if len(genomes_annotated) > 1:
            return [recommendation(
                CHECK, "alignment", "",
                f"the scheme annotates {len(genomes_annotated)} reference genomes "
                f"({', '.join(genomes_annotated)}), so --genome must say which one the alignment "
                f"is against; the genomes were not checked")]

    try:
        mapping = AlignmentMap(next(iter(alignment.values())), reference)
    except ReferenceError as error:
        return [problem(CHECK, "alignment", "",
                        f"the alignment's first record is not the reference genome: {error}")]

    names = list(alignment)
    clusters, unmatched, ambiguous = resolve_clusters(clusters, names)
    placement = Placement(scheme, clusters, names)
    placement.missing_genomes = sorted(set(placement.missing_genomes) | set(unmatched))
    if ambiguous:
        findings.append(problem(
            CHECK, "ambiguous-genome", "",
            f"{len(ambiguous)} genomes of the genome map match more than one record of the "
            f"alignment, and were left out: {_names(ambiguous)}"))

    if placement.unknown_clusters:
        findings.append(problem(
            CHECK, "unknown-cluster", "",
            f"the genome map places {len(placement.unplaced)} genomes at "
            f"{len(placement.unknown_clusters)} nodes the scheme does not hold: "
            f"{_names(placement.unknown_clusters)}. Where they belong in this scheme is not "
            f"known, so they were left out of every check, as genomes neither below nor "
            f"outside a node. If export --prune left those nodes out, the genomes file it "
            f"writes places their genomes at the nearest node kept"))
    if placement.missing_genomes:
        findings.append(recommendation(
            CHECK, "missing-genome", "",
            f"the genome map names {len(placement.missing_genomes)} genomes the alignment does "
            f"not hold, and they were left out: {_names(placement.missing_genomes)}"))
    if placement.ungrouped:
        findings.append(recommendation(
            CHECK, "ungrouped", "",
            f"{len(placement.ungrouped)} genomes of the alignment are placed at no node, and "
            f"count as outside every node: {_names(placement.ungrouped)}"))

    snp_findings, failed = check_scheme_snps(scheme, alignment, mapping, placement, genome)
    findings += snp_findings
    findings += check_genome_paths(scheme, alignment, mapping, placement, genome,
                                   genome_conflicts)

    if new_clades:
        canonical, clades, by_node, by_genome, sites = scan_alignment(alignment, mapping,
                                                                      placement)
        findings += check_extension(scheme, mapping, placement, canonical, clades, by_node,
                                    by_genome, failed, sites, min_clade_size, min_clade_snps)
        logging.info(
            f"Read {sites} variable sites over {len(names)} genomes, of which "
            f"{len(placement.node_of)} are placed in the scheme")
    else:
        logging.info(
            f"Checked the scheme's SNPs in {len(names)} genomes, of which "
            f"{len(placement.node_of)} are placed in the scheme; clades the scheme does not "
            f"resolve were not looked for")
    return findings
