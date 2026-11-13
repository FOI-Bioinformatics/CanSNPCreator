"""How a scheme types the genomes of an alignment, and what that says about it.

The audit's other checks against an alignment start from a placement: a genome
map or a tree says where each genome belongs, and the SNPs are judged against
it. This module needs no placement. It types every genome of the alignment
with the scheme's own canSNPs, the way a typing tool would, and measures the
scheme by the result, so that a scheme can be evaluated on genomes it was not
made from.

Each genome is typed under two rules:

``single``
    The rule of CanSNPer2 2.2.0, which reads one canSNP per node, the one of
    rank 1. The call is the deepest node whose whole chain of nodes from the
    root is derived, when that chain holds at least ``min_required_hits``
    nodes. A node without a base at its canSNP breaks the chain.
``set``
    A node is derived when more of its canSNPs carry the derived base than
    the ancestral one, ancestral when fewer, and unknown otherwise. The call
    is the deepest derived node with no ancestral node above it; a node of
    unknown state may lie on the way.

With one canSNP per node the two rules differ only in how they treat missing
bases and short chains. With several, the set rule also outvotes a single
discordant base. The difference between the two columns of the report is
therefore what a set of canSNPs per node would gain on these genomes.

Three tables and a summary are written:

- per genome: the call under each rule and the bases behind it;
- per canSNP: how often it is missing, and how often it disagrees with the
  other canSNPs of its node or with the path of the genome. The evaluate
  subcommand reads this table to rank candidates by how they behaved
  (``--observed``);
- per node: how many genomes pass through it and how much support the set
  rule had there.

Genomes that contradict the scheme at ``genome_conflicts`` nodes or more are
typed and listed, but are left out of the per canSNP disagreement counts: on
1,229 Brucella genomes half of all such contradictions came from 5% of the
genomes, which says more about those assemblies than about the SNPs.
"""

import logging
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

from cansnpscorer.helpers.data_utils import MINIMUM_SNPS_PER_NODE
from cansnpscorer.helpers.names import resolve_clusters
from cansnpscorer.model import AlignmentMap, Reference
from cansnpscorer.scheme import Scheme, SchemeSnp

BASES = frozenset("ACGT")

DERIVED, ANCESTRAL, MISSING, OTHER = "D", "A", "N", "O"

GENOME_COLUMNS = ("genome", "single_call", "single_depth", "single_stop", "set_call",
                  "set_depth", "derived_on_path", "missing_on_path", "ancestral_on_path",
                  "derived_off_path", "discordant", "expected", "single_outcome",
                  "set_outcome")
SNP_COLUMNS = ("node", "rank", "genome", "position", "derived_base", "ancestral_base",
               "genomes", "called", "missing", "derived", "ancestral", "other",
               "discordant_in_set", "ancestral_on_path", "derived_off_path")
NODE_COLUMNS = ("node", "depth", "snps", "single_calls", "set_calls", "through",
                "all_derived", "unknown_on_path", "split", "min_margin", "single_stops")
SUMMARY_COLUMNS = ("metric", "single", "set")

#: How a call compares with the node a genome map expects.
SAME, ABOVE, BELOW, ELSEWHERE, NO_CALL = "same", "above", "below", "other", "none"


class TypingError(ValueError):
    """The scheme cannot be used to type this alignment."""


@dataclass
class GenomeTyping:
    """One genome typed under both rules."""

    genome: str
    single_call: str | None = None
    single_stop: str = ""
    """Why the single rule ended above the set call: the node its chain broke
    at and that node's state, or ``min-hits`` for a complete chain too short."""
    set_call: str | None = None
    derived_on_path: int = 0
    missing_on_path: int = 0
    ancestral_on_path: int = 0
    derived_off_path: int = 0
    discordant: bool = False
    expected: str | None = None
    states: dict[tuple[str, int], str] = field(default_factory=dict, repr=False)
    """The state of each canSNP in this genome, by (node, index within node)."""

    @property
    def conflicts(self) -> int:
        return self.ancestral_on_path + self.derived_off_path


@dataclass
class TypingReport:
    """Every genome typed, and the scheme measured by it."""

    genomes: list[GenomeTyping]
    snp_rows: list[list[str]]
    node_rows: list[list[str]]
    summary: list[tuple[str, str, str]]
    depth: dict[str, int]
    lineage: dict[str, list[str]] = field(default_factory=dict, repr=False)

    def genome_rows(self) -> list[list[str]]:
        rows = []
        for typed in self.genomes:
            rows.append([
                typed.genome,
                typed.single_call or "NA",
                str(self.depth.get(typed.single_call, 0)) if typed.single_call else "0",
                typed.single_stop,
                typed.set_call or "NA",
                str(self.depth.get(typed.set_call, 0)) if typed.set_call else "0",
                str(typed.derived_on_path), str(typed.missing_on_path),
                str(typed.ancestral_on_path), str(typed.derived_off_path),
                "yes" if typed.discordant else "no",
                typed.expected or "",
                outcome(typed.single_call, typed.expected, self.lineage)
                if typed.expected else "",
                outcome(typed.set_call, typed.expected, self.lineage)
                if typed.expected else "",
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


def single_call(first: Mapping[str, str], lineage: Mapping[str, Sequence[str]],
                min_required_hits: int) -> str | None:
    """The call of CanSNPer2 2.2.0 from the state of each node's first canSNP."""
    complete = [node for node, state in first.items()
                if state == DERIVED and len(lineage[node]) >= min_required_hits
                and all(first.get(above) == DERIVED for above in lineage[node])]
    return deepest(complete, lineage)


def set_call(states: Mapping[str, str], lineage: Mapping[str, Sequence[str]]) -> str | None:
    """The deepest derived node with no ancestral node above it."""
    reachable = [node for node, state in states.items()
                 if state == DERIVED
                 and not any(states.get(above) == ANCESTRAL for above in lineage[node])]
    return deepest(reachable, lineage)


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
                   genome_conflicts: int = 3) -> TypingReport:
    """Type every genome of the alignment with the scheme, under both rules.

    Args:
        scheme: The scheme as read from its source files.
        alignment: The whole genome alignment, reference first.
        reference: The reference genome the scheme is annotated against.
        genome: Which of the scheme's reference genomes the alignment is
            against, when it names more than one.
        clusters: Where a genome map expects the genomes, when one is known;
            each call is then compared with it.
        min_required_hits: The nodes the single rule requires on a chain.
        genome_conflicts: A genome contradicting the scheme at this many
            nodes or more is marked discordant.

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

    expected: dict[str, str] = {}
    if clusters:
        resolved, _unmatched, _ambiguous = resolve_clusters(clusters, list(alignment))
        for node, names in resolved.items():
            if node in hierarchy.parents:
                for name in names:
                    expected.setdefault(name, node)

    typed: list[GenomeTyping] = []
    for name, sequence in alignment.items():
        result = GenomeTyping(genome=name, expected=expected.get(name))
        per_node: dict[str, list[str]] = {}
        for node, snps in by_node.items():
            states = [MISSING if column is None else snp_state(sequence[column], snp)
                      for snp, column in snps]
            per_node[node] = states
            for index, state in enumerate(states):
                result.states[(node, index)] = state

        nodes = {node: node_state(states) for node, states in per_node.items()}
        first = {node: (states[0] if states[0] in (DERIVED, ANCESTRAL) else MISSING)
                 for node, states in per_node.items()}
        result.set_call = set_call(nodes, lineage)
        result.single_call = single_call(first, lineage, min_required_hits)

        path = set(lineage[result.set_call]) if result.set_call else set()
        for node, state in nodes.items():
            if node in path:
                if state == DERIVED:
                    result.derived_on_path += 1
                elif state == ANCESTRAL:  # cannot arise: such a node blocks the path
                    result.ancestral_on_path += 1
                else:
                    result.missing_on_path += 1
            elif state == DERIVED and result.set_call:
                result.derived_off_path += 1
        # Nodes of the path that carry no canSNP at all are missing too.
        result.missing_on_path += sum(1 for node in path if node not in nodes)
        result.discordant = result.conflicts >= genome_conflicts
        result.single_stop = _single_stop(result, first, lineage, min_required_hits)
        typed.append(result)

    report = TypingReport(genomes=typed, snp_rows=[], node_rows=[], summary=[], depth=depth,
                          lineage=lineage)
    report.snp_rows = _snp_rows(by_node, typed, lineage)
    report.node_rows = _node_rows(hierarchy.nodes, by_node, typed, lineage)
    report.summary = _summary(typed, depth, lineage, by_node)
    return report


def _single_stop(result: GenomeTyping, first: Mapping[str, str],
                 lineage: Mapping[str, Sequence[str]], min_required_hits: int) -> str:
    """Why the single rule did not reach the set call, when it did not."""
    if result.set_call is None or result.single_call == result.set_call:
        return ""
    if result.single_call and result.set_call in lineage[result.single_call]:
        return ""  # the first canSNPs lead further than the majority of the sets
    for node in reversed(lineage[result.set_call]):  # from the root down
        state = first.get(node, MISSING)
        if state != DERIVED:
            return f"{node}:{'ancestral' if state == ANCESTRAL else 'missing'}"
    if len(lineage[result.set_call]) < min_required_hits:
        return "min-hits"
    return "ambiguous"


def _snp_rows(by_node, typed: Sequence[GenomeTyping],
              lineage: Mapping[str, Sequence[str]]) -> list[list[str]]:
    rows = []
    for node, snps in by_node.items():
        # How many of the node's canSNPs are derived and ancestral in each
        # genome, counted once, so that a node of thousands of candidates
        # costs no more per canSNP than a node of one.
        totals = []
        for result in typed:
            states = [result.states[(node, index)] for index in range(len(snps))]
            totals.append((states.count(DERIVED), states.count(ANCESTRAL)))
        for index, (snp, _column) in enumerate(snps):
            counts: Counter = Counter()
            for result, (derived, ancestral) in zip(typed, totals, strict=True):
                state = result.states[(node, index)]
                counts[state] += 1
                if result.discordant or state not in (DERIVED, ANCESTRAL):
                    continue
                # The other canSNPs of the node: its totals without this one.
                if (state == DERIVED and ancestral > derived - 1) or (state == ANCESTRAL and derived > ancestral - 1):
                    counts["discordant"] += 1
                if result.set_call:
                    on_path = node in lineage[result.set_call]
                    if on_path and state == ANCESTRAL:
                        counts["ancestral_on_path"] += 1
                    elif not on_path and state == DERIVED:
                        counts["derived_off_path"] += 1
            called = counts[DERIVED] + counts[ANCESTRAL] + counts[OTHER]
            rows.append([
                node, str(snp.rank), snp.genome, str(snp.position), snp.derived, snp.ancestral,
                str(len(typed)), str(called), str(counts[MISSING]), str(counts[DERIVED]),
                str(counts[ANCESTRAL]), str(counts[OTHER]), str(counts["discordant"]),
                str(counts["ancestral_on_path"]), str(counts["derived_off_path"]),
            ])
    return rows


def _node_rows(nodes: Sequence[str], by_node, typed: Sequence[GenomeTyping],
               lineage: Mapping[str, Sequence[str]]) -> list[list[str]]:
    rows = []
    single_calls = Counter(result.single_call for result in typed)
    set_calls = Counter(result.set_call for result in typed)
    stops = Counter(result.single_stop.partition(":")[0] for result in typed
                    if ":" in result.single_stop)
    for node in nodes:
        snps = by_node.get(node, [])
        through = all_derived = unknown = split = 0
        margin: int | None = None
        for result in typed:
            states = [result.states[(node, index)] for index in range(len(snps))]
            derived, ancestral = states.count(DERIVED), states.count(ANCESTRAL)
            if derived and ancestral:
                split += 1
            if not result.set_call or node not in lineage[result.set_call]:
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
            "" if margin is None else str(margin), str(stops[node]),
        ])
    return rows


def _summary(typed: Sequence[GenomeTyping], depth: Mapping[str, int],
             lineage: Mapping[str, Sequence[str]], by_node) -> list[tuple[str, str, str]]:
    total = len(typed)
    singles = [result.single_call for result in typed]
    sets = [result.set_call for result in typed]

    def mean_depth(calls: Sequence[str | None]) -> str:
        return f"{sum(depth[call] for call in calls if call) / total:.2f}" if total else "0"

    sizes = sorted(len(snps) for snps in by_node.values())
    rows = [
        ("genomes", str(total), str(total)),
        ("canSNPs per node, fewest to most", f"{sizes[0]}-{sizes[-1]}", f"{sizes[0]}-{sizes[-1]}"),
        ("typed", str(sum(1 for call in singles if call)), str(sum(1 for call in sets if call))),
        ("not typed", str(sum(1 for call in singles if not call)),
         str(sum(1 for call in sets if not call))),
        ("mean depth of call", mean_depth(singles), mean_depth(sets)),
        ("nodes called at", str(len({call for call in singles if call})),
         str(len({call for call in sets if call}))),
        ("single call above the set call", str(sum(
            1 for result in typed if result.set_call and result.single_call != result.set_call
            and (result.single_call is None
                 or result.single_call in lineage[result.set_call]))), ""),
        ("genomes with a derived node off their path", "",
         str(sum(1 for result in typed if result.derived_off_path))),
        ("discordant genomes", "", str(sum(1 for result in typed if result.discordant))),
    ]
    with_expected = [result for result in typed if result.expected]
    if with_expected:
        rows.append(("genomes with an expected node", str(len(with_expected)),
                     str(len(with_expected))))
        for kind in (SAME, ABOVE, BELOW, ELSEWHERE, NO_CALL):
            rows.append((
                f"call {kind} the expected node" if kind in (ABOVE, BELOW)
                else {SAME: "call is the expected node", ELSEWHERE: "call on another branch",
                      NO_CALL: "no call"}[kind],
                str(sum(1 for r in with_expected
                        if outcome(r.single_call, r.expected, lineage) == kind)),
                str(sum(1 for r in with_expected
                        if outcome(r.set_call, r.expected, lineage) == kind))))
    return rows


def write_table(path: str, header: Sequence[str], rows: Sequence[Sequence[str]]) -> None:
    with open(path, "w", newline="\n") as handle:
        handle.write("\t".join(header) + "\n")
        handle.writelines("\t".join(row) + "\n" for row in rows)


def write_report(report: TypingReport, prefix: str) -> list[str]:
    """Write the four tables as PREFIX_genomes.tsv, _snps.tsv, _nodes.tsv and
    _summary.tsv. Returns the paths written."""
    tables = (("_genomes.tsv", GENOME_COLUMNS, report.genome_rows()),
              ("_snps.tsv", SNP_COLUMNS, report.snp_rows),
              ("_nodes.tsv", NODE_COLUMNS, report.node_rows),
              ("_summary.tsv", SUMMARY_COLUMNS, report.summary))
    written = []
    for suffix, header, rows in tables:
        write_table(prefix + suffix, header, rows)
        written.append(prefix + suffix)
    for metric, single, by_set in report.summary:
        logging.info(f"Typing: {metric}: " + ", ".join(
            f"{rule} {value}" for rule, value in (("single", single), ("set", by_set)) if value))
    logging.info(f"Wrote the typing report to {', '.join(written)}")
    return written


@dataclass
class Observed:
    """How one canSNP behaved in the genomes of a typing report."""

    genomes: int
    missing: int
    discordant: int
    """Genomes where it disagreed with the other canSNPs of its node or with
    the genome's path."""


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
        lacking = [column for column in SNP_COLUMNS if column not in header]
        if lacking:
            raise ValueError(
                f"{path}: has no column {', '.join(lacking)}; it should be the _snps.tsv table "
                f"of audit --typing-report")
        at = {column: header.index(column) for column in SNP_COLUMNS}
        for number, line in enumerate(handle, start=2):
            fields = line.rstrip("\n").split("\t")
            if not line.strip():
                continue
            try:
                key = (fields[at["node"]], int(fields[at["position"]]))
                genomes.add(fields[at["genome"]])
                observed[key] = Observed(
                    genomes=max(int(fields[at["genomes"]]), 1),
                    missing=int(fields[at["missing"]]) + int(fields[at["other"]]),
                    discordant=max(int(fields[at["discordant_in_set"]]),
                                   int(fields[at["ancestral_on_path"]])
                                   + int(fields[at["derived_off_path"]])))
            except (IndexError, ValueError):
                raise ValueError(f"{path}: line {number} is not a row of the table") from None
    if not observed:
        raise ValueError(f"{path}: holds no canSNPs")
    if len(genomes) > 1:
        raise ValueError(
            f"{path}: holds canSNPs of {len(genomes)} reference genomes "
            f"({', '.join(sorted(genomes))}); run the typing report with --genome")
    return observed
