"""How a scheme types the genomes of an alignment, and what that says about it.

The audit's other checks against an alignment start from a placement: a genome
map or a tree says where each genome belongs, and the SNPs are judged against
it. This module needs no placement. It types every genome of the alignment
with the scheme's own canSNPs, the way a typing tool would, and measures the
scheme by the result, so that a scheme can be evaluated on genomes it was not
made from.

Each genome is typed under two rules:

``single``
    The rule of CanSNPer2 2.2.1, which reads one canSNP per node, the one of
    rank 1. The call is the deepest derived node whose chain of nodes from
    the root holds no ancestral node and no node of another base, holds at
    least ``min_required_hits`` derived nodes that are at least
    ``strictness`` of the chain, and whose parent is derived. A node without
    a base at its canSNP is a gap in the chain, not a break.
``set``
    A node is derived when more of its canSNPs carry the derived base than
    the ancestral one, ancestral when fewer, and unknown otherwise. A derived
    node with no ancestral node above it can be the call; a node of unknown
    state may lie on the way to it. Of such nodes the call is the one with
    the most derived nodes on its path from the root, and the deepest of
    those. Counting the support, not the depth alone, keeps one homoplastic
    canSNP below a run of unknown nodes from outweighing a well supported
    path. Two best supported nodes on different branches give no call.

With one canSNP per node the two rules differ only in how they treat missing
bases and short chains. With several, the set rule also outvotes a single
discordant base. The difference between the two columns of the report is
therefore what a set of canSNPs per node would gain on these genomes.

Four tables, a summary and a record of the run are written:

- per genome: the call under each rule and the bases behind it;
- per canSNP: how often it is missing, and how often it disagrees with the
  other canSNPs of its node or with the path of the genome. The evaluate
  subcommand reads this table to rank candidates by how they behaved
  (``--observed``);
- per node: how many genomes pass through it and how much support the set
  rule had there;
- per rank: the same counts summed over the canSNPs of each rank, which
  shows whether the first canSNP of a node behaves better than the others.

The reference genome and the outgroups are typed and listed like any genome,
but the summary counts neither: the reference has a base at every canSNP by
construction, and an outgroup is expected to get no call. A field that does
not apply holds ``NA`` in every table.

Genomes that contradict the scheme at ``genome_conflicts`` nodes or more are
typed and listed, but are left out of the per canSNP disagreement counts: on
1,229 Brucella genomes half of all such contradictions came from 5% of the
genomes, which says more about those assemblies than about the SNPs.
"""

import logging
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace

from cansnpscorer.helpers.data_utils import MINIMUM_SNPS_PER_NODE
from cansnpscorer.helpers.names import resolve_clusters, resolve_names
from cansnpscorer.model import AlignmentMap, Reference
from cansnpscorer.scheme import Scheme, SchemeSnp

BASES = frozenset("ACGT")

DERIVED, ANCESTRAL, MISSING, OTHER = "D", "A", "N", "O"

GENOME_COLUMNS = ("genome", "role", "single_call", "single_depth", "single_stop_node",
                  "single_stop_reason", "set_call",
                  "set_depth", "derived_on_path", "missing_on_path", "ancestral_on_path",
                  "derived_off_path", "discordant", "expected", "single_outcome",
                  "set_outcome")
SNP_COLUMNS = ("node", "rank", "reference_genome", "position", "contig", "contig_position",
               "derived_base", "ancestral_base",
               "genomes", "called", "missing", "derived", "ancestral", "other",
               "discordant_in_set", "ancestral_on_path", "derived_off_path", "derived_set",
               "set_apart")
#: The columns evaluate --observed needs. Tables written before 0.4 name the
#: reference genome ``genome`` and hold no contig; they are still read.
OBSERVED_COLUMNS = ("node", "position", "derived_base", "ancestral_base", "missing", "derived",
                    "ancestral", "other", "discordant_in_set", "ancestral_on_path",
                    "derived_off_path")
RANK_COLUMNS = ("rank", "snps", "genomes", "called", "missing", "missing_pct",
                "discordant_in_set", "ancestral_on_path", "derived_off_path",
                "snps_disagreeing", "snps_disagreeing_pct")
NODE_COLUMNS = ("node", "depth", "snps", "single_calls", "set_calls", "through",
                "all_derived", "unknown_on_path", "split", "min_margin", "single_stops")
SUMMARY_COLUMNS = ("metric", "single", "single_pct", "set", "set_pct")
RUN_COLUMNS = ("item", "value")

NA = "NA"
#: A node of more canSNPs than this has only its first ones labelled by the
#: genomes they are derived in, which bounds the pairwise comparison.
MAX_SET_LABELLED = 200
#: What a record of the alignment is to the scheme.
REFERENCE, OUTGROUP, ISOLATE = "reference", "outgroup", "isolate"

#: How a call compares with the node a genome map expects.
SAME, ABOVE, BELOW, ELSEWHERE, NO_CALL = "same", "above", "below", "other", "none"


class TypingError(ValueError):
    """The scheme cannot be used to type this alignment."""


@dataclass
class GenomeTyping:
    """One genome typed under both rules."""

    genome: str
    role: str = ISOLATE
    single_call: str | None = None
    single_stop: tuple[str, str] = ("", "")
    """Why the single rule did not give the set call, as (node, reason): the
    node its chain broke at with ``ancestral`` or ``other-base``; the parent
    of the set call without a base with ``missing`` or ``other-reference``;
    ``min-hits`` or ``strictness`` for too few derived nodes; ``other-branch``
    with the node it called instead; ``ambiguous`` for two confirmed chains
    of equal depth."""
    set_call: str | None = None
    path: frozenset[str] = frozenset()
    """The nodes from the root to the set call; where two branches tie and no
    call is made, the nodes they share."""
    derived_on_path: int = 0
    missing_on_path: int = 0
    ancestral_on_path: int = 0
    derived_off_path: int = 0
    discordant: bool = False
    expected: str | None = None
    states: str = field(default="", repr=False)
    """The state of each canSNP in this genome, one character each, in the
    order of the report's canSNPs. A string and not a table by canSNP: with
    20,000 candidates and 1,200 genomes a table took about 2 GB."""

    @property
    def conflicts(self) -> int:
        return self.ancestral_on_path + self.derived_off_path


@dataclass
class TypingReport:
    """Every genome typed, and the scheme measured by it."""

    genomes: list[GenomeTyping]
    snp_rows: list[list[str]]
    node_rows: list[list[str]]
    summary: list[tuple[str, str, str, str, str]]
    depth: dict[str, int]
    lineage: dict[str, list[str]] = field(default_factory=dict, repr=False)
    rank_rows: list[list[str]] = field(default_factory=list)
    stops: list[tuple[str, str, int]] = field(default_factory=list)
    """Where the single rule's chains broke, as (node, reason, genomes), the
    most frequent first."""
    unreached: list[str] = field(default_factory=list)
    """Nodes carrying canSNPs that no counted genome passes through, so that
    these genomes say nothing about them."""

    def genome_rows(self) -> list[list[str]]:
        rows = []
        for typed in self.genomes:
            rows.append([
                typed.genome, typed.role,
                typed.single_call or NA,
                str(self.depth[typed.single_call]) if typed.single_call else NA,
                *(part or NA for part in typed.single_stop),
                typed.set_call or NA,
                str(self.depth[typed.set_call]) if typed.set_call else NA,
                str(typed.derived_on_path), str(typed.missing_on_path),
                str(typed.ancestral_on_path), str(typed.derived_off_path),
                "yes" if typed.discordant else "no",
                typed.expected or NA,
                outcome(typed.single_call, typed.expected, self.lineage)
                if typed.expected else NA,
                outcome(typed.set_call, typed.expected, self.lineage)
                if typed.expected else NA,
            ])
        return rows


def outcome(call: str | None, expected: str, lineage: Mapping[str, Sequence[str]]) -> str:
    """How a call compares with the node expected for the genome."""
    if call is None:
        return NO_CALL
    if call == expected:
        return SAME
    if call in lineage.get(expected, ()):
        return ABOVE
    if expected in lineage.get(call, ()):
        return BELOW
    return ELSEWHERE


def _lineages(hierarchy) -> dict[str, list[str]]:
    """Each node with every node above it, itself first."""
    lineage: dict[str, list[str]] = {}
    for node in hierarchy.nodes:  # parents come before their children
        parent = hierarchy.parents[node]
        lineage[node] = [node, *(lineage[parent] if parent is not None else [])]
    return lineage


def snp_state(base: str, snp: SchemeSnp) -> str:
    base = base.upper()
    if base == snp.derived.upper():
        return DERIVED
    if base == snp.ancestral.upper():
        return ANCESTRAL
    return OTHER if base in BASES else MISSING


def node_state(states: Sequence[str]) -> str:
    """The state of a node from those of its canSNPs: the majority of the
    derived and ancestral ones, and missing when they are equal or absent."""
    derived, ancestral = states.count(DERIVED), states.count(ANCESTRAL)
    if derived > ancestral:
        return DERIVED
    if ancestral > derived:
        return ANCESTRAL
    return MISSING


def deepest(candidates: Sequence[str], lineage: Mapping[str, Sequence[str]]) -> str | None:
    """The deepest of some nodes, or None when two share the greatest depth:
    they then lie on different branches, and neither is the call."""
    if not candidates:
        return None
    depth = max(len(lineage[node]) for node in candidates)
    lowest = [node for node in candidates if len(lineage[node]) == depth]
    return lowest[0] if len(lowest) == 1 else None


#: CanSNPer2's default --strictness: the fraction of the canSNPs on a chain
#: that must be derived.
STRICTNESS = 0.7


def confirm_chain(node: str, first: Mapping[str, str], lineage: Mapping[str, Sequence[str]],
                  min_required_hits: int, strictness: float = STRICTNESS) -> tuple[str, str]:
    """Whether CanSNPer2 2.2.1 confirms the chain from the root to a node.

    Each node of the chain is derived, ancestral, of another base, or
    without a base. An ancestral node or one of another base breaks the
    chain; a node without a base is a gap in it. The chain is confirmed when
    it is not broken, at least ``min_required_hits`` of its nodes are
    derived, the derived fraction of its nodes is at least ``strictness``,
    and the node directly above the last one has a base (``NewickTree.
    _confirm_path`` of CanSNPer2 2.2.1).

    Returns:
        ("", "") when confirmed; otherwise the reason and the node it
        concerns: (node, "ancestral"), (node, "other-base") for the topmost
        break, (parent, "missing") for a parent without a base, or ("",
        "min-hits") and ("", "strictness") for too few derived nodes.
    """
    chain = lineage[node]
    for above in reversed(chain):  # from the root down
        state = first.get(above, MISSING)
        if state == ANCESTRAL:
            return above, "ancestral"
        if state == OTHER:
            return above, "other-base"
    derived = sum(1 for above in chain if first.get(above) == DERIVED)
    if len(chain) > 1 and first.get(chain[1]) != DERIVED:
        return chain[1], "missing"
    if derived < min_required_hits:
        return "", "min-hits"
    if derived < strictness * len(chain):
        return "", "strictness"
    return "", ""


def single_call(first: Mapping[str, str], lineage: Mapping[str, Sequence[str]],
                min_required_hits: int, strictness: float = STRICTNESS) -> str | None:
    """The call of CanSNPer2 2.2.1 from the state of each node's first canSNP:
    the deepest derived node whose chain ``confirm_chain`` confirms, when the
    genome has at least ``min_required_hits`` derived nodes. Two confirmed
    nodes at the greatest depth give no call."""
    derived = [node for node, state in first.items() if state == DERIVED]
    if len(derived) < min_required_hits:
        return None
    confirmed = [node for node in derived
                 if confirm_chain(node, first, lineage, min_required_hits, strictness) == ("", "")]
    return deepest(confirmed, lineage)


def set_path(states: Mapping[str, str], lineage: Mapping[str, Sequence[str]]
             ) -> tuple[str | None, frozenset[str]]:
    """The call of the set rule and the nodes from the root to it.

    A derived node with no ancestral node above it can be the call. The one
    with the most derived nodes on its lineage is taken, and the deepest of
    those. When the best supported nodes lie on different branches there is
    no call, and the path is the part of their lineages they share, so that
    what the genome carries beyond it still counts as contradicting the
    scheme.
    """
    reachable = [node for node, state in states.items()
                 if state == DERIVED
                 and not any(states.get(above) == ANCESTRAL for above in lineage[node])]
    if not reachable:
        return None, frozenset()
    support = {node: sum(1 for above in lineage[node] if states.get(above) == DERIVED)
               for node in reachable}
    most = max(support.values())
    best = [node for node in reachable if support[node] == most]
    lowest = max(best, key=lambda node: len(lineage[node]))
    if all(node in lineage[lowest] for node in best):
        return lowest, frozenset(lineage[lowest])
    shared = set(lineage[best[0]])
    for node in best[1:]:
        shared &= set(lineage[node])
    return None, frozenset(shared)


def set_call(states: Mapping[str, str], lineage: Mapping[str, Sequence[str]]) -> str | None:
    """The call of the set rule alone; see :func:`set_path`."""
    return set_path(states, lineage)[0]


def derived_off(states: Mapping[str, str], lineage: Mapping[str, Sequence[str]],
                call: str | None, path: frozenset[str]) -> int:
    """How many derived nodes contradict the path of a genome.

    With a call these are the derived nodes off its path. Where branches tie
    and there is no call, one of the tied branches may still be the genome's
    own, so the count is the fewest nodes that contradict any one of them:
    the derived nodes beyond the shared path, less those of the tied branch
    that holds most of them.
    """
    off = [node for node, state in states.items() if state == DERIVED and node not in path]
    if call is not None or not off or not path:
        return len(off)
    beyond = set(off)
    own = max(sum(1 for above in lineage[node] if above in beyond) for node in off)
    return len(off) - own


def scheme_columns(scheme: Scheme, mapping: AlignmentMap, genome: str | None
                   ) -> dict[str, list[tuple[SchemeSnp, int | None]]]:
    """The canSNPs of each node, best rank first, with their alignment column."""
    by_node: dict[str, list[tuple[SchemeSnp, int | None]]] = {}
    for snp in sorted(scheme.checked_snps, key=lambda snp: (snp.rank, snp.line)):
        if genome is not None and snp.genome != genome:
            continue
        if snp.node not in scheme.hierarchy.parents:
            continue
        column = None
        if 1 <= snp.position <= mapping.reference.length:
            column = mapping.to_column(snp.position - 1)
        by_node.setdefault(snp.node, []).append((snp, column))
    return by_node


def type_alignment(scheme: Scheme, alignment: Mapping[str, str], reference: Reference,
                   genome: str | None = None,
                   clusters: Mapping[str, Sequence[str]] | None = None,
                   min_required_hits: int = MINIMUM_SNPS_PER_NODE,
                   genome_conflicts: int = 3,
                   strictness: float = STRICTNESS,
                   outgroups: Sequence[str] = ()) -> TypingReport:
    """Type every genome of the alignment with the scheme, under both rules.

    Args:
        scheme: The scheme as read from its source files.
        alignment: The whole genome alignment, reference first.
        reference: The reference genome the scheme is annotated against.
        genome: Which of the scheme's reference genomes the alignment is
            against, when it names more than one.
        clusters: Where a genome map expects the genomes, when one is known;
            each call is then compared with it.
        min_required_hits: The derived nodes the single rule requires on a
            chain.
        genome_conflicts: A genome contradicting the scheme at this many
            nodes or more is marked discordant.
        outgroups: Genomes known to lie outside the scheme, by the names a
            genome map would give them. They are typed, but not counted.
        strictness: The fraction of a chain's nodes the single rule requires
            to be derived, as CanSNPer2's --strictness.

    Raises:
        TypingError: If the scheme has no single tree, or does not say which
            reference genome the alignment is against.
        ReferenceError: If the alignment's first record is not the reference.
    """
    if scheme.hierarchy is None:
        raise TypingError("the tree file is not a single rooted tree, so no genome can be typed")
    annotated = sorted({snp.genome for snp in scheme.snps})
    if genome is None and len(annotated) > 1:
        raise TypingError(
            f"the scheme annotates {len(annotated)} reference genomes ({', '.join(annotated)}), "
            f"so --genome must say which one the alignment is against")

    mapping = AlignmentMap(next(iter(alignment.values())), reference)
    hierarchy = scheme.hierarchy
    lineage = _lineages(hierarchy)
    depth = {node: len(above) for node, above in lineage.items()}
    by_node = scheme_columns(scheme, mapping, genome)
    if not by_node:
        raise TypingError("the scheme holds no canSNP at a node of its tree for this genome")
    # The single rule reads the canSNP of rank 1 and no other.
    has_first = {node: snps[0][0].rank == 1 for node, snps in by_node.items()}
    # Nodes whose canSNPs are all annotated on another reference genome cannot
    # be read in this alignment. CanSNPer2 aligns to every reference, so the
    # single rule's column is not its call below such a node.
    elsewhere = {snp.node for snp in scheme.checked_snps} - set(by_node)
    if elsewhere:
        logging.warning(
            f"{len(elsewhere)} nodes carry canSNPs on other reference genomes only, and count "
            f"as without a base here: {', '.join(sorted(elsewhere)[:5])}. Below them the "
            f"single rule's call is not what CanSNPer2 gives, which aligns to every reference")

    expected: dict[str, str] = {}
    if clusters:
        resolved, _unmatched, _ambiguous = resolve_clusters(clusters, list(alignment))
        for node, names in resolved.items():
            if node in hierarchy.parents:
                for name in names:
                    expected.setdefault(name, node)

    outside: set[str] = set()
    if outgroups:
        matched, _unmatched, _ambiguous = resolve_names(list(outgroups), list(alignment))
        outside = set(matched.values())
    first_record = next(iter(alignment))

    typed: list[GenomeTyping] = []
    for name, sequence in alignment.items():
        role = REFERENCE if name == first_record else OUTGROUP if name in outside else ISOLATE
        result = GenomeTyping(genome=name, role=role, expected=expected.get(name))
        per_node: dict[str, list[str]] = {}
        for node, snps in by_node.items():
            per_node[node] = [MISSING if column is None else snp_state(sequence[column], snp)
                              for snp, column in snps]
        result.states = "".join("".join(states) for states in per_node.values())

        nodes = {node: node_state(states) for node, states in per_node.items()}
        # CanSNPer2 2.2.1 tells a third base, which breaks a chain, from no
        # base, which is a gap in it.
        first = {node: (states[0] if has_first[node] else MISSING)
                 for node, states in per_node.items()}
        result.set_call, result.path = set_path(nodes, lineage)
        result.single_call = single_call(first, lineage, min_required_hits, strictness)

        path = result.path
        for node, state in nodes.items():
            if node in path:
                if state == DERIVED:
                    result.derived_on_path += 1
                elif state == ANCESTRAL:  # cannot arise: such a node blocks the path
                    result.ancestral_on_path += 1
                else:
                    result.missing_on_path += 1
        result.derived_off_path = derived_off(nodes, lineage, result.set_call, path)
        # Nodes of the path that carry no canSNP at all are missing too.
        result.missing_on_path += sum(1 for node in path if node not in nodes)
        result.discordant = result.conflicts >= genome_conflicts
        result.single_stop = _single_stop(result, first, lineage, min_required_hits,
                                          elsewhere, strictness)
        typed.append(result)

    report = TypingReport(genomes=typed, snp_rows=[], node_rows=[], summary=[], depth=depth,
                          lineage=lineage)
    # Where the canSNPs of each node begin in a genome's string of states.
    starts: dict[str, int] = {}
    total = 0
    for node, snps in by_node.items():
        starts[node] = total
        total += len(snps)
    counted = [result for result in typed if result.role == ISOLATE]
    report.snp_rows = _snp_rows(by_node, starts, typed, reference)
    report.node_rows = _node_rows(hierarchy.nodes, by_node, starts, typed, lineage)
    report.rank_rows = _rank_rows(report.snp_rows)
    report.stops = [(node, reason, count) for (node, reason), count in sorted(
        Counter(result.single_stop for result in counted
                if result.single_stop[1] in BROKEN).items(),
        key=lambda item: (-item[1], item[0]))]
    on_a_path = {node for result in counted for node in result.path}
    report.unreached = [node for node in hierarchy.nodes
                        if node in by_node and node not in on_a_path]
    report.summary = _summary(typed, depth, lineage, by_node, report)
    return report


#: The reasons that name the node where the single rule's chain broke.
BROKEN = ("missing", "ancestral", "other-base", "other-reference")


def _single_stop(result: GenomeTyping, first: Mapping[str, str],
                 lineage: Mapping[str, Sequence[str]], min_required_hits: int,
                 elsewhere: frozenset[str] | set[str] = frozenset(),
                 strictness: float = STRICTNESS) -> tuple[str, str]:
    """Why the single rule did not give the set call, as (node, reason).

    The reason is why ``confirm_chain`` does not confirm the chain to the
    set call; a parent without a base because its canSNP is annotated on
    another reference genome is "other-reference". When the chain is
    confirmed the reason is "min-hits" if the genome has fewer derived nodes
    than required, and otherwise "ambiguous": another confirmed node lies as
    deep or deeper.
    """
    if result.set_call is None or result.single_call == result.set_call:
        return "", ""
    if result.single_call and result.set_call in lineage[result.single_call]:
        return "", ""  # the first canSNPs lead further than the majority of the sets
    if result.single_call and result.single_call not in lineage[result.set_call]:
        return result.single_call, "other-branch"
    node, reason = confirm_chain(result.set_call, first, lineage, min_required_hits,
                                 strictness)
    if reason == "missing" and node in elsewhere:
        return node, "other-reference"
    if reason:
        return node, reason
    if sum(1 for state in first.values() if state == DERIVED) < min_required_hits:
        return "", "min-hits"
    return "", "ambiguous"


def derived_sets(start: int, count: int, typed: Sequence[GenomeTyping]) -> list[int]:
    """Label the canSNPs of one node by the genomes they are derived in.

    Two canSNPs share a label when the genomes derived in one and ancestral
    in the other number no more than ``width_tolerance`` of the genomes, the
    tolerance of the count rule; a genome without a base may be in either
    state. Labels are given greedily from the canSNP derived in the most
    genomes, each canSNP taking the first label whose first canSNP it agrees
    with, so that label 1 is the set of most genomes. A canSNP without a base
    in more genomes than the tolerance, which ``width_groups`` leaves out,
    founds no label others may join: it takes the first label it agrees
    with, or else a label of its own. Genomes marked discordant are left
    out, as in the disagreement counts. A canSNP derived in none of the
    remaining genomes gets 0, as do the canSNPs beyond the first
    ``MAX_SET_LABELLED`` of a node; 0 is written as an empty field. Missing
    bases here are counted over the non-discordant genomes, whereas
    ``width_groups`` counts them over all genomes.
    """
    judged = [result.states for result in typed if not result.discordant]
    tolerance = width_tolerance(len(judged))
    width = min(count, MAX_SET_LABELLED)
    sets = []
    for index in range(start, start + width):
        column = [states[index] for states in judged]
        sets.append((frozenset(g for g, state in enumerate(column) if state == DERIVED),
                     frozenset(g for g, state in enumerate(column) if state == ANCESTRAL)))
    order = sorted((index for index in range(width) if sets[index][0]),
                   key=lambda index: (-len(sets[index][0]), index))
    covered, sparse = [], []
    for index in order:
        missing = len(judged) - len(sets[index][0]) - len(sets[index][1])
        (covered if missing <= tolerance else sparse).append(index)
    labels = [0] * count
    first: list[tuple[frozenset[int], frozenset[int]]] = []

    def agreeing(index: int) -> int:
        derived, ancestral = sets[index]
        for label, (other_derived, other_ancestral) in enumerate(first, start=1):
            if (len(derived & other_ancestral) + len(ancestral & other_derived)
                    <= tolerance):
                return label
        return 0

    for index in covered:
        labels[index] = agreeing(index)
        if not labels[index]:
            first.append(sets[index])
            labels[index] = len(first)
    last = len(first)
    for index in sparse:
        label = agreeing(index)
        if not label:
            last += 1
            label = last
        labels[index] = label
    return labels


def set_apart(start: int, count: int, typed: Sequence[GenomeTyping],
              labels: Sequence[int]) -> list[int | None]:
    """How far each canSNP of one node lies from the narrowest of its sets.

    The sets are those of ``derived_sets``, as ``width_groups`` reads them
    back: a canSNP without a base in more genomes than the tolerance is in
    none, and those derived in no genome form a set of their own, the
    narrowest. The narrowest set that holds at least two canSNPs gives each
    genome a state, the one most of its canSNPs have there, and none where
    they are as often derived as ancestral. A canSNP lies as far from it as
    the genomes in which it has the other state; a genome where either has
    no base counts for neither. Genomes marked discordant are left out, as in
    ``derived_sets``.

    Returns:
        One distance per canSNP, or None for every canSNP when the node has
        no set of two, or holds more canSNPs than ``MAX_SET_LABELLED``.
    """
    if count > MAX_SET_LABELLED:
        return [None] * count
    tolerance = width_tolerance(len(typed))
    members: dict[int, list[int]] = {}
    for index in range(count):
        column = [result.states[start + index] for result in typed]
        derived = column.count(DERIVED)
        if derived + column.count(ANCESTRAL) == 0:
            continue
        if len(column) - derived - column.count(ANCESTRAL) > tolerance:
            continue
        members.setdefault(labels[index] if derived else 0, []).append(index)
    sets = [label for label in sorted(members, key=lambda label: (not label, label))
            if len(members[label]) > 1]
    if not sets:
        return [None] * count
    narrowest = members[sets[-1]]
    judged = [result.states[start:start + count] for result in typed if not result.discordant]
    reference = []
    for states in judged:
        votes = Counter(states[index] for index in narrowest)
        reference.append(DERIVED if votes[DERIVED] > votes[ANCESTRAL]
                         else ANCESTRAL if votes[ANCESTRAL] > votes[DERIVED] else MISSING)
    apart = []
    for index in range(count):
        apart.append(sum(1 for states, state in zip(judged, reference, strict=True)
                         if state != MISSING and states[index] in (DERIVED, ANCESTRAL)
                         and states[index] != state))
    return apart


def _snp_rows(by_node, starts: Mapping[str, int], typed: Sequence[GenomeTyping],
              reference: Reference) -> list[list[str]]:
    rows = []
    crowded = [node for node, snps in by_node.items() if len(snps) > MAX_SET_LABELLED]
    if crowded:
        logging.info(
            f"Typing: {len(crowded)} nodes hold more than {MAX_SET_LABELLED} canSNPs; only "
            f"their first {MAX_SET_LABELLED} are labelled in derived_set, and evaluate "
            f"--observed groups their canSNPs by count: {', '.join(crowded[:5])}"
            + (" and others" if len(crowded) > 5 else ""))
    for node, snps in by_node.items():
        start = starts[node]
        labels = derived_sets(start, len(snps), typed)
        distances = set_apart(start, len(snps), typed, labels)
        # How many of the node's canSNPs are derived and ancestral in each
        # genome, counted once, so that a node of thousands of candidates
        # costs no more per canSNP than a node of one.
        totals = []
        for result in typed:
            states = result.states[start:start + len(snps)]
            totals.append((states.count(DERIVED), states.count(ANCESTRAL)))
        for index, (snp, _column) in enumerate(snps):
            counts: Counter = Counter()
            for result, (derived, ancestral) in zip(typed, totals, strict=True):
                state = result.states[start + index]
                counts[state] += 1
                if result.discordant or state not in (DERIVED, ANCESTRAL):
                    continue
                # The other canSNPs of the node: its totals without this one.
                if (state == DERIVED and ancestral > derived - 1) or (
                        state == ANCESTRAL and derived > ancestral - 1):
                    counts["discordant"] += 1
                on_path = node in result.path
                if on_path and state == ANCESTRAL:
                    counts["ancestral_on_path"] += 1
                elif not on_path and state == DERIVED:
                    counts["derived_off_path"] += 1
            called = counts[DERIVED] + counts[ANCESTRAL] + counts[OTHER]
            contig, local = NA, NA
            if 1 <= snp.position <= reference.length:
                found, offset = reference.locate(snp.position - 1)
                contig, local = found.name, str(offset + 1)
            rows.append([
                node, str(snp.rank), snp.genome, str(snp.position), contig, local,
                snp.derived, snp.ancestral,
                str(len(typed)), str(called), str(counts[MISSING]), str(counts[DERIVED]),
                str(counts[ANCESTRAL]), str(counts[OTHER]), str(counts["discordant"]),
                str(counts["ancestral_on_path"]), str(counts["derived_off_path"]),
                str(labels[index]) if labels[index] else "",
                "" if distances[index] is None else str(distances[index]),
            ])
    return rows


def percent(part: int, whole: int) -> str:
    return f"{100 * part / whole:.1f}" if whole else NA


def _rank_rows(snp_rows: Sequence[Sequence[str]]) -> list[list[str]]:
    """The counts of the canSNP table summed over the canSNPs of each rank."""
    at = {column: SNP_COLUMNS.index(column) for column in SNP_COLUMNS}
    summed = ("genomes", "called", "missing", "discordant_in_set", "ancestral_on_path",
              "derived_off_path")
    by_rank: dict[int, Counter] = {}
    for row in snp_rows:
        counts = by_rank.setdefault(int(row[at["rank"]]), Counter())
        counts["snps"] += 1
        for column in summed:
            counts[column] += int(row[at[column]])
        if any(int(row[at[column]]) for column in summed[3:]):
            counts["disagreeing"] += 1
    rows = []
    for rank, counts in sorted(by_rank.items()):
        rows.append([
            str(rank), str(counts["snps"]), str(counts["genomes"]), str(counts["called"]),
            str(counts["missing"]), percent(counts["missing"], counts["genomes"]),
            str(counts["discordant_in_set"]), str(counts["ancestral_on_path"]),
            str(counts["derived_off_path"]), str(counts["disagreeing"]),
            percent(counts["disagreeing"], counts["snps"]),
        ])
    return rows


def _node_rows(nodes: Sequence[str], by_node, starts: Mapping[str, int],
               typed: Sequence[GenomeTyping],
               lineage: Mapping[str, Sequence[str]]) -> list[list[str]]:
    rows = []
    single_calls = Counter(result.single_call for result in typed)
    set_calls = Counter(result.set_call for result in typed)
    stops = Counter(result.single_stop[0] for result in typed
                    if result.single_stop[1] in BROKEN)
    for node in nodes:
        snps = by_node.get(node, [])
        start = starts.get(node, 0)
        through = all_derived = unknown = split = 0
        margin: int | None = None
        for result in typed:
            states = result.states[start:start + len(snps)]
            derived, ancestral = states.count(DERIVED), states.count(ANCESTRAL)
            if derived and ancestral:
                split += 1
            if node not in result.path:
                continue
            through += 1
            if derived > ancestral:
                margin = derived - ancestral if margin is None else min(margin,
                                                                       derived - ancestral)
                if derived == len(snps):
                    all_derived += 1
            else:
                unknown += 1
        rows.append([
            node, str(len(lineage[node])), str(len(snps)), str(single_calls[node]),
            str(set_calls[node]), str(through), str(all_derived), str(unknown), str(split),
            NA if margin is None else str(margin), str(stops[node]),
        ])
    return rows


#: How many nodes the summary names where the single rule stopped.
NAMED_STOPS = 5


def _summary(everyone: Sequence[GenomeTyping], depth: Mapping[str, int],
             lineage: Mapping[str, Sequence[str]], by_node,
             report: "TypingReport") -> list[tuple[str, str, str, str, str]]:
    typed = [result for result in everyone if result.role == ISOLATE]
    total = len(typed)
    singles = [result.single_call for result in typed]
    sets = [result.set_call for result in typed]

    def mean_depth(calls: Sequence[str | None]) -> str:
        made = [depth[call] for call in calls if call]
        return f"{sum(made) / len(made):.2f}" if made else NA

    def compared(result: GenomeTyping) -> str:
        """The single call beside the set call, where the set rule made one."""
        single, by_set = result.single_call, result.set_call
        if single == by_set:
            return "same"
        if single is None:
            return "none"
        if single in lineage[by_set]:
            return "above"
        return "below" if by_set in lineage[single] else "other"

    def text(metric: str, single: str = NA, by_set: str = NA):
        return (metric, single, NA, by_set, NA)

    def count(metric: str, single: int | None, by_set: int | None, whole: int = total):
        return (metric,
                NA if single is None else str(single),
                NA if single is None else percent(single, whole),
                NA if by_set is None else str(by_set),
                NA if by_set is None else percent(by_set, whole))

    beside = Counter(compared(result) for result in typed if result.set_call)
    only_single = sum(1 for result in typed if result.single_call and not result.set_call)
    reasons = Counter(result.single_stop[1] for result in typed if result.single_stop[1])
    references = [result for result in everyone if result.role == REFERENCE]
    outgroups = [result for result in everyone if result.role == OUTGROUP]

    sizes = sorted(len(snps) for snps in by_node.values())
    rows = [
        text("records of the alignment", str(len(everyone)), str(len(everyone))),
        text("reference genome, typed at",
             *((references[0].single_call or NA, references[0].set_call or NA)
               if references else (NA, NA))),
        text("outgroups, not counted below", str(len(outgroups)), str(len(outgroups))),
        text("outgroups given a call",
             str(sum(1 for result in outgroups if result.single_call)),
             str(sum(1 for result in outgroups if result.set_call))),
        text("genomes", str(total), str(total)),
        text("nodes with canSNPs", str(len(by_node)), str(len(by_node))),
        text("canSNPs per node, fewest to most", f"{sizes[0]}-{sizes[-1]}",
             f"{sizes[0]}-{sizes[-1]}"),
        count("typed", sum(1 for call in singles if call), sum(1 for call in sets if call)),
        count("not typed", sum(1 for call in singles if not call),
              sum(1 for call in sets if not call)),
        text("mean depth of the calls made", mean_depth(singles), mean_depth(sets)),
        text("nodes called at", str(len({call for call in singles if call})),
             str(len({call for call in sets if call}))),
        text("nodes no genome passes through", NA, str(len(report.unreached))),
        count("single call is the set call", beside["same"], None),
        count("single call above the set call", beside["above"], None),
        count("single call below the set call", beside["below"], None),
        count("single call on another branch than the set call", beside["other"], None),
        count("no single call where the set rule calls", beside["none"], None),
        count("single call where the set rule gives none", only_single, None),
        count("single rule stopped at a canSNP without a base", reasons["missing"], None),
        count("single rule stopped at an ancestral canSNP", reasons["ancestral"], None),
        count("single rule stopped at a node on another reference",
              reasons["other-reference"], None),
        count("single rule stopped for too few nodes on the chain", reasons["min-hits"], None),
        count("single rule stopped at two chains of equal depth", reasons["ambiguous"], None),
        text("nodes where the single rule stopped most often",
             "; ".join(f"{node} ({found} {reason})"
                       for node, reason, found in report.stops[:NAMED_STOPS]) or NA),
        count("genomes with a derived node off their path", None,
              sum(1 for result in typed if result.derived_off_path)),
        count("discordant genomes", None, sum(1 for result in typed if result.discordant)),
    ]
    with_expected = [result for result in typed if result.expected]
    if with_expected:
        rows.append(text("genomes with an expected node", str(len(with_expected)),
                         str(len(with_expected))))
        for kind in (SAME, ABOVE, BELOW, ELSEWHERE, NO_CALL):
            rows.append(count(
                f"call {kind} the expected node" if kind in (ABOVE, BELOW)
                else {SAME: "call is the expected node", ELSEWHERE: "call on another branch",
                      NO_CALL: "no call"}[kind],
                sum(1 for r in with_expected
                    if outcome(r.single_call, r.expected, lineage) == kind),
                sum(1 for r in with_expected
                    if outcome(r.set_call, r.expected, lineage) == kind),
                len(with_expected)))
    return rows


def write_table(path: str, header: Sequence[str], rows: Sequence[Sequence[str]]) -> None:
    with open(path, "w", newline="\n") as handle:
        handle.write("\t".join(header) + "\n")
        handle.writelines("\t".join(row) + "\n" for row in rows)


def write_report(report: TypingReport, prefix: str,
                 run: Sequence[tuple[str, str]] = ()) -> list[str]:
    """Write the tables as PREFIX_genomes.tsv, _snps.tsv, _nodes.tsv,
    _ranks.tsv and _summary.tsv, and with ``run``, what the report was made
    from, as PREFIX_run.tsv. Returns the paths written."""
    tables = [("_genomes.tsv", GENOME_COLUMNS, report.genome_rows()),
              ("_snps.tsv", SNP_COLUMNS, report.snp_rows),
              ("_nodes.tsv", NODE_COLUMNS, report.node_rows),
              ("_ranks.tsv", RANK_COLUMNS, report.rank_rows),
              ("_summary.tsv", SUMMARY_COLUMNS, report.summary)]
    if run:
        tables.append(("_run.tsv", RUN_COLUMNS, [list(item) for item in run]))
    written = []
    for suffix, header, rows in tables:
        write_table(prefix + suffix, header, rows)
        written.append(prefix + suffix)
    for metric, single, single_pct, by_set, set_pct in report.summary:
        logging.info(f"Typing: {metric}: " + ", ".join(
            f"{rule} {value}" + (f" ({share}%)" if share != NA else "")
            for rule, value, share in (("single", single, single_pct), ("set", by_set, set_pct))
            if value != NA))
    if report.unreached:
        logging.info(
            f"Typing: {len(report.unreached)} nodes lie on the path of no genome, so these "
            f"genomes do not test their canSNPs: {', '.join(report.unreached[:5])}"
            + (" and others" if len(report.unreached) > 5 else ""))
    logging.info(f"Wrote the typing report to {', '.join(written)}")
    return written


@dataclass
class Observed:
    """How one canSNP behaved in the genomes of a typing report."""

    genomes: int
    """Genomes in which it had the derived or the ancestral base; 0 when it
    was called in none, which is no observation at all."""
    missing: int
    discordant: int
    """Genomes where it disagreed with the other canSNPs of its node or with
    the genome's path."""
    derived: str = ""
    ancestral: str = ""
    derived_genomes: int = 0
    """Genomes in which it had the derived base."""
    derived_set: int = 0
    """Which set of genomes it is derived in among the canSNPs of its node,
    1 for the set of most genomes; 0 when unknown, as in tables written
    before the column was added."""
    set_apart: int | None = None
    """Genomes in which it has the other state than the narrowest set of two
    or more canSNPs of its node; None when unknown, as in tables written
    before the column was added."""


def read_observed(path: str) -> dict[tuple[str, int], Observed]:
    """Read PREFIX_snps.tsv, keyed by (node, 1-based position over the
    concatenated reference).

    Raises:
        ValueError: If the table is not one, or holds canSNPs of more than
            one reference genome, whose positions cannot be told apart.
    """
    observed: dict[tuple[str, int], Observed] = {}
    genomes: set[str] = set()
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        lacking = [column for column in OBSERVED_COLUMNS if column not in header]
        named = "reference_genome" if "reference_genome" in header else "genome"
        if named not in header:
            lacking.append("reference_genome")
        if lacking:
            raise ValueError(
                f"{path}: has no column {', '.join(lacking)}; it should be the _snps.tsv table "
                f"of audit --typing-report")
        at = {column: header.index(column) for column in (*OBSERVED_COLUMNS, named)}
        # Tables of 0.4.0 and earlier have no derived_set column.
        labelled = header.index("derived_set") if "derived_set" in header else None
        distant = header.index("set_apart") if "set_apart" in header else None
        for number, line in enumerate(handle, start=2):
            fields = line.rstrip("\n").split("\t")
            if not line.strip():
                continue
            try:
                key = (fields[at["node"]], int(fields[at["position"]]))
                genomes.add(fields[at[named]])
                row = Observed(
                    genomes=int(fields[at["derived"]]) + int(fields[at["ancestral"]]),
                    missing=int(fields[at["missing"]]) + int(fields[at["other"]]),
                    discordant=max(int(fields[at["discordant_in_set"]]),
                                   int(fields[at["ancestral_on_path"]])
                                   + int(fields[at["derived_off_path"]])),
                    derived=fields[at["derived_base"]].upper(),
                    ancestral=fields[at["ancestral_base"]].upper(),
                    derived_genomes=int(fields[at["derived"]]),
                    derived_set=(int(fields[labelled]) if labelled is not None
                                 and labelled < len(fields) and fields[labelled] else 0),
                    set_apart=(int(fields[distant]) if distant is not None
                               and distant < len(fields) and fields[distant] else None))
                earlier = observed.get(key)
                if earlier is not None:
                    # Listed twice, as a node with two ranks at one position:
                    # the worse of the two is kept, not whichever came last.
                    row = Observed(min(row.genomes, earlier.genomes),
                                   max(row.missing, earlier.missing),
                                   max(row.discordant, earlier.discordant),
                                   row.derived, row.ancestral,
                                   min(row.derived_genomes, earlier.derived_genomes),
                                   row.derived_set if row.derived_set == earlier.derived_set
                                   else 0,
                                   # set_apart counted both rows as members of
                                   # their set; read as one, the sets may differ.
                                   None)
                observed[key] = row
            except (IndexError, ValueError):
                raise ValueError(f"{path}: line {number} is not a row of the table") from None
    if not observed:
        raise ValueError(f"{path}: holds no canSNPs")
    if len(genomes) > 1:
        raise ValueError(
            f"{path}: holds canSNPs of {len(genomes)} reference genomes "
            f"({', '.join(sorted(genomes))}); run the typing report with --genome")
    return observed


def width_tolerance(genomes: int) -> int:
    """How far the derived counts of two canSNPs may lie apart and still count
    as the same width: one genome in a hundred, and at least two."""
    return max(2, genomes // 100)


def width_groups(rows: Sequence[Observed]) -> list[list[Observed]]:
    """The groups of ``_width_groups``, without how they were told apart."""
    return _width_groups(rows)[0]


def _width_groups(rows: Sequence[Observed]) -> tuple[list[list[Observed]], bool]:
    """The canSNPs of one node grouped by the genomes they are derived in,
    widest group first.

    When every canSNP derived in some genome carries a ``derived_set``
    label, the canSNPs of one label are one group, and those derived in no
    genome a last group. Otherwise, as for a table without the column, the
    groups are told apart by the number of genomes alone: two canSNPs are of
    one group when their counts are joined by steps of at most the
    tolerance. In both cases a canSNP that lacks a base in more genomes than
    the tolerance has no known width and is in no group.

    Returns:
        The groups, and whether they were told apart by the labels.
    """
    groups: list[list[Observed]] = []
    known = [row for row in rows
             if row.genomes and row.missing <= width_tolerance(row.genomes + row.missing)]
    if any(row.derived_set for row in known) and all(
            row.derived_set for row in known if row.derived_genomes):
        by_label: dict[int, list[Observed]] = {}
        for row in sorted(known, key=lambda row: -row.derived_genomes):
            by_label.setdefault(row.derived_set if row.derived_genomes else 0, []).append(row)
        return [by_label[label] for label in sorted(by_label, key=lambda label: (not label,
                                                                                   label))], True
    for row in sorted(known, key=lambda row: -row.derived_genomes):
        tolerance = width_tolerance(row.genomes + row.missing)
        if groups and groups[-1][-1].derived_genomes - row.derived_genomes <= tolerance:
            groups[-1].append(row)
        else:
            groups.append([row])
    return groups, False


def judge_by_narrowest_group(
        observed: Mapping[tuple[str, int], Observed]
) -> tuple[dict[tuple[str, int], Observed], list[str]]:
    """Judge the canSNPs of nodes that hold no majority by their narrowest group.

    The typing report counts a canSNP as disagreeing where it differs from
    the majority of its node's canSNPs. That presumes the canSNPs of a node
    mark one group of genomes. They do among the genomes the scheme was made
    from, but genomes typed later may join the node's branch partway along
    it and carry its earlier mutations only, so that the canSNPs fall into
    nested groups of different width. When one group holds more than half of
    the canSNPs the majority is that group, and the count stands. When none
    does, the majority follows no group, and the canSNPs it favours are a
    mixture.

    For such a node the canSNPs are judged against one group instead: the
    narrowest that holds at least two canSNPs, whose genomes carry every
    mutation of the branch and so are the clade of the node's own genomes.

    Groups are told apart by the genomes their canSNPs are derived in, as
    the ``derived_set`` column of the table labels them, so that two groups
    of equal size on different genomes are two groups. The disagreement of
    each canSNP then becomes the number of genomes in which it has the other
    state than that group, as the ``set_apart`` column gives it.

    A table without the ``derived_set`` column, or a node where a canSNP
    derived in some genome lacks a label, tells the groups apart by their
    size alone, which nested groups always differ in. Where the groups are
    told apart by labels but a canSNP lacks a distance, as one listed twice
    at a position does, the groups stand and only the distance falls back.
    Measured by counts, the disagreement of each canSNP is the number of
    genomes by which its derived count must differ from the group's most
    common count. In both cases a genome without a base allows for either
    state.

    Returns:
        The observations, changed for those nodes, and the nodes.
    """
    by_node: dict[str, list[tuple[tuple[str, int], Observed]]] = {}
    for key, row in observed.items():
        by_node.setdefault(key[0], []).append((key, row))
    judged = dict(observed)
    changed = []
    for node, found in by_node.items():
        groups, by_set = _width_groups([row for _key, row in found])
        if len(groups) < 2 or max(len(group) for group in groups) * 2 > len(found):
            continue
        narrowest = next((group for group in reversed(groups) if len(group) > 1), None)
        if narrowest is None:
            continue
        counts = Counter(row.derived_genomes for row in narrowest)
        width = min(counts, key=lambda count: (-counts[count], count))
        by_set = by_set and all(row.set_apart is not None for _key, row in found
                                if row.genomes)
        for key, row in found:
            if not row.genomes:
                continue
            if by_set:
                apart = row.set_apart
            else:
                apart = max(0, row.derived_genomes - width,
                            width - row.derived_genomes - row.missing)
            judged[key] = replace(row, discordant=apart)
        changed.append(node)
    return judged, changed
