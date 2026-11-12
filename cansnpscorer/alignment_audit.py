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
from typing import Dict, FrozenSet, List, Mapping, Optional, Sequence, Set, Tuple

from cansnpscorer.findings import CHECK, EXTEND, MODIFY, Finding, problem, recommendation
from cansnpscorer.helpers.data_utils import snp_extractor
from cansnpscorer.helpers.names import resolve_clusters
from cansnpscorer.model import AlignmentMap, Reference, ReferenceError
from cansnpscorer.scheme import Scheme

BASES = frozenset("ACGT")

#: At most this many genome names are written into one finding.
SHOWN = 5
#: At most this many new clades are reported one by one; the rest are counted.
REPORTED_CLADES = 50


def _names(genomes: Sequence[str]) -> str:
    shown = sorted(genomes)[:SHOWN]
    more = len(genomes) - len(shown)
    return ", ".join(shown) + (f" and {more} more" if more > 0 else "")


def _where(reference: Reference, column_position: int) -> str:
    contig, local = reference.locate(column_position)
    return f"{contig.name}:{local + 1}"


class Placement:
    """Where the genome map puts each genome, and the isolates below each node."""

    def __init__(self, scheme: Scheme, clusters: Mapping[str, Sequence[str]],
                 genomes: Sequence[str]) -> None:
        hierarchy = scheme.hierarchy
        self.hierarchy = hierarchy
        self.genomes: List[str] = list(genomes)
        known = set(genomes)

        self.unknown_clusters = sorted(set(clusters) - set(hierarchy.parents))
        self.missing_genomes = sorted(
            {name for node, names in clusters.items() for name in names} - known)

        usable = {
            node: [name for name in names if name in known]
            for node, names in clusters.items() if node in hierarchy
        }
        self.node_of: Dict[str, str] = {}
        for node, names in usable.items():
            for name in names:
                self.node_of.setdefault(name, node)

        below = hierarchy.isolates(usable)
        self.below: Dict[str, FrozenSet[str]] = {
            node: frozenset(names) for node, names in below.items()}
        self.ungrouped = sorted(known - set(self.node_of))
        self.node_by_set: Dict[FrozenSet[str], str] = {}
        # Parents come before children, so the deepest node of a set is kept;
        # a node with one child holds the same isolates as that child.
        for node in hierarchy.nodes:
            if self.below[node]:
                self.node_by_set[self.below[node]] = node

    def lowest_holding(self, carriers: FrozenSet[str]) -> Optional[str]:
        """The deepest node whose isolates include all of these genomes."""
        first = next(iter(carriers))
        node: Optional[str] = self.node_of.get(first)
        while node is not None and not carriers <= self.below[node]:
            node = self.hierarchy.parents[node]
        return node

    def fits(self, carriers: FrozenSet[str], node: str) -> bool:
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
              genomes: Sequence[str]) -> Dict[str, List[str]]:
    by_base: Dict[str, List[str]] = {}
    for genome in genomes:
        base = alignment[genome][column].upper()
        by_base.setdefault(base if base in BASES else "missing", []).append(genome)
    return by_base


def check_scheme_snps(scheme: Scheme, alignment: Mapping[str, str], mapping: AlignmentMap,
                      placement: Placement, genome: Optional[str]) -> Tuple[List[Finding],
                                                                           Set[str]]:
    """Whether each SNP of the scheme is canonical in these genomes.

    Returns the findings and the nodes whose SNP is not canonical.
    """
    findings: List[Finding] = []
    failed: Set[str] = set()
    untested: List[str] = []
    canonical = 0

    for snp in scheme.snps:
        if genome is not None and snp.genome != genome:
            continue
        if snp.node not in placement.below or snp.position < 1:
            continue
        if snp.position > mapping.reference.length:
            continue  # reported by the reference check
        inside = placement.below[snp.node]
        if not inside:
            untested.append(snp.node)
            continue

        column = mapping.to_column(snp.position - 1)
        if column is None:
            findings.append(recommendation(
                CHECK, "not-aligned", snp.node,
                f"no column of the alignment covers position {snp.position} "
                f"({_where(mapping.reference, snp.position - 1)}), so the SNP cannot be checked "
                f"in these genomes"))
            continue

        outside = [name for name in placement.genomes if name not in inside]
        if not outside:
            # Every genome sits below this node, which is the root in practice.
            # A SNP all of them carry separates nothing in this alignment, and
            # its polarity cannot be told without a genome outside the node.
            findings.append(recommendation(
                CHECK, "no-outgroup", snp.node,
                f"every genome of the alignment sits below this node, so nothing tests whether "
                f"its SNP at {snp.position} marks it, nor which base is ancestral. Add an "
                f"outgroup genome to the alignment, placed at no node"))
            continue

        within = _bases_at(alignment, column, sorted(inside))
        beyond = _bases_at(alignment, column, outside)
        derived, ancestral = snp.derived.upper(), snp.ancestral.upper()
        where = f"at {snp.position} ({_where(mapping.reference, snp.position - 1)})"

        called_inside = {base for base in within if base != "missing"}
        not_derived = [name for base, names in within.items()
                       if base not in (derived, "missing") for name in names]
        derived_outside = beyond.get(derived, [])
        missing = within.get("missing", []) + beyond.get("missing", [])

        if called_inside == {ancestral} and derived_outside:
            failed.add(snp.node)
            findings.append(problem(
                MODIFY, "polarity", snp.node,
                f"the isolates below the node carry '{ancestral}' {where}, which the scheme gives "
                f"as the ancestral base, and {len(derived_outside)} genomes outside it carry the "
                f"derived '{derived}'. The two columns appear swapped"))
            continue

        if not_derived:
            failed.add(snp.node)
            findings.append(problem(
                MODIFY, "not-canonical", snp.node,
                f"{len(not_derived)} of {len(inside)} isolates below the node do not carry the "
                f"derived '{derived}' {where}: {_names(not_derived)}. CanSNPer2 would place "
                f"them above this node"))

        if derived_outside:
            failed.add(snp.node)
            placed = Counter(placement.node_of.get(name, "no node") for name in derived_outside)
            nodes = ", ".join(f"{node} ({count})" for node, count in placed.most_common(SHOWN))
            findings.append(problem(
                MODIFY, "homoplasy", snp.node,
                f"{len(derived_outside)} genomes outside the node carry the derived "
                f"'{derived}' {where}: {_names(derived_outside)}; they are placed at {nodes}. "
                f"The SNP does not mark this branch alone"))

        if missing:
            findings.append(recommendation(
                CHECK, "missing-data", snp.node,
                f"{len(missing)} genomes have no base {where}: {_names(missing)}. The SNP cannot "
                f"be called in them, which is how an accessory or poorly aligned region shows "
                f"itself"))

        if snp.node not in failed:
            canonical += 1

    if untested:
        findings.append(recommendation(
            CHECK, "untested", ", ".join(sorted(untested)[:SHOWN]),
            f"{len(untested)} nodes have no genome placed at or below them, so their SNPs were "
            f"not checked"))

    logging.info(
        f"{canonical} SNPs of the scheme are canonical in these genomes, and {len(failed)} "
        f"are not")
    return findings, failed


def scan_alignment(alignment: Mapping[str, str], mapping: AlignmentMap,
                   placement: Placement) -> Tuple[Dict[str, List[Tuple[int, str, str]]],
                                                  Dict[Tuple[str, FrozenSet[str]],
                                                       List[Tuple[int, str, str]]],
                                                  Counter, Counter, int]:
    """Sort every variable site by what it says about the scheme.

    Only sites where every genome has one of two bases are used, as the
    extract subcommand uses them.

    Returns:
        - canonical sites per node, as (position, derived, ancestral)
        - sites per candidate new clade, keyed by (parent node, genomes)
        - conflicting sites per node
        - conflicting sites per genome
        - the number of sites read
    """
    canonical: Dict[str, List[Tuple[int, str, str]]] = {}
    clades: Dict[Tuple[str, FrozenSet[str]], List[Tuple[int, str, str]]] = {}
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
        carriers: Dict[str, Set[str]] = {}
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
                    canonical: Dict[str, List[Tuple[int, str, str]]],
                    clades: Dict[Tuple[str, FrozenSet[str]], List[Tuple[int, str, str]]],
                    conflicts_by_node: Counter, conflicts_by_genome: Counter, failed: Set[str],
                    sites: int, min_clade_size: int = 2, min_clade_snps: int = 1
                    ) -> List[Finding]:
    """Replacements for failed SNPs, nodes these genomes do not support, and new clades."""
    findings: List[Finding] = []
    reference = mapping.reference

    def described(entries: Sequence[Tuple[int, str, str]]) -> str:
        shown = [f"{position + 1} ({_where(reference, position)}, {ancestral}>{derived})"
                 for position, derived, ancestral in entries[:3]]
        return "; ".join(shown) + (" and others" if len(entries) > 3 else "")

    for node in scheme.hierarchy.nodes:
        if not placement.below.get(node):
            continue
        found = canonical.get(node, [])
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
            findings.append(recommendation(
                CHECK, "single-site", node,
                "one site of the alignment separates exactly the isolates below this node, the "
                "scheme's own. A sequencing error in that one base moves a sample off the "
                "branch, and nothing else in these genomes supports it"))

    ranked = sorted(clades.items(), key=lambda item: (-len(item[1]), item[0][0],
                                                      sorted(item[0][1])))

    # Each grouping fits the scheme on its own, but two of them can still
    # contradict each other: genomes shared without one holding the other. No
    # tree holds both, so the better supported one is kept and the other is
    # counted with the conflicting sites, as homoplasy or noise.
    accepted: List[FrozenSet[str]] = []
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
                            genome: Optional[str] = None, min_clade_size: int = 2,
                            min_clade_snps: int = 1) -> List[Finding]:
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
    """
    findings: List[Finding] = []

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
            f"the genome map places genomes at {len(placement.unknown_clusters)} nodes the scheme "
            f"does not hold, and they were left out: {_names(placement.unknown_clusters)}"))
    if placement.missing_genomes:
        findings.append(problem(
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

    canonical, clades, by_node, by_genome, sites = scan_alignment(alignment, mapping, placement)
    findings += check_extension(scheme, mapping, placement, canonical, clades, by_node,
                                by_genome, failed, sites, min_clade_size, min_clade_snps)

    logging.info(
        f"Read {sites} variable sites over {len(names)} genomes, of which "
        f"{len(names) - len(placement.ungrouped)} are placed in the scheme")
    return findings
