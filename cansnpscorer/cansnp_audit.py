"""Audit of an existing CanSNPer2-data scheme.

The validation script of the CanSNPer2-data repository answers whether a scheme
is well formed: whether the three files parse, whether every node carries a SNP
and whether the accessions look right. It does not answer the question asked of
a published scheme in practice, which is whether the scheme is any good: which
branches cannot be called at all, which SNPs can never be called derived, which
positions look as though they were taken from a different genome version, and
where a better SNP is available.

This module answers that question. It reads a scheme back in, optionally beside
the reference genome it is annotated against and the scores from the evaluate
subcommand, and reports findings of two kinds:

- a problem, which makes part of the scheme wrong or impossible to type, and
- a recommendation, which leaves the scheme usable but says how it would be
  improved.

Each finding carries the action it asks for: ``modify`` an existing SNP or
node, ``extend`` the scheme with SNPs it does not have, or ``check`` something
this tool cannot settle on its own.

Two of the checks come from schemes already in the repository. The ported
Coxiella scheme had its positions taken from an earlier version of the
reference assembly, which the reference base check finds directly; the ported
Brucella scheme has branches two nodes deep, which CanSNPer2 cannot call with
its default of three required hits.
"""

import itertools
import logging
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace

from cansnpscorer import __version__, vcf
from cansnpscorer.alignment_audit import GENOME_CONFLICTS, check_against_alignment
from cansnpscorer.cansnp_evaluator import read_masks
from cansnpscorer.cansnp_export import GENBANK, REFSEQ
from cansnpscorer.findings import (
    CHECK,
    EXTEND,
    MODIFY,
    PROBLEM,
    RECOMMENDATION,
    Finding,
    problem,
    recommendation,
)
from cansnpscorer.helpers.data_utils import MINIMUM_SNPS_PER_NODE
from cansnpscorer.helpers.file_utils import (
    ColumnRow,
    read_alignment,
    read_alignment_columns,
    read_contig_aliases,
    read_genome_map,
    reference_parser,
)
from cansnpscorer.helpers.scoring_methods import ScoredSnp
from cansnpscorer.hierarchy import HierarchyError
from cansnpscorer.model import AlignmentMap, Reference, ReferenceError
from cansnpscorer.report import write_markdown
from cansnpscorer.scheme import Scheme, SchemeError, read_scheme
from cansnpscorer.tree_tools import place, read_genome_tree, read_members
from cansnpscorer.typing_report import (
    NA,
    STRICTNESS,
    TypingError,
    TypingReport,
    scheme_columns,
    type_alignment,
    write_report,
)

BASES = frozenset("ACGT")


def check_structure(scheme: Scheme) -> list[Finding]:
    """The tree, and whether every node carries exactly one SNP."""
    findings: list[Finding] = []

    for node, parents in scheme.conflicts:
        findings.append(
            problem(
                MODIFY,
                "two-parents",
                node,
                f"the node is placed under more than one parent ({', '.join(parents)}). The "
                f"CanSNPer2 tree reader does not refuse this: it creates a second node of the same "
                f"name and attaches the SNP to one of them",
            )
        )

    by_node = scheme.snps_by_node
    nodes = scheme.nodes

    for node in nodes:
        if ";" in node:
            findings.append(
                problem(
                    MODIFY,
                    "semicolon",
                    node,
                    "the node name holds a semicolon. The CanSNPer2 tree reader turns tabs into "
                    "semicolons before splitting, so the name is read as two columns",
                )
            )

    for node in nodes:
        found = by_node.get(node, [])
        if not found:
            findings.append(
                problem(
                    EXTEND,
                    "no-snp",
                    node,
                    "the node is in the tree but has no SNP, so it can never be called. Every node "
                    "of a CanSNPer2 scheme carries exactly one",
                )
            )
        elif len(found) > 1:
            findings.append(
                problem(
                    MODIFY,
                    "repeated-snp",
                    node,
                    f"{len(found)} SNPs annotate the node (lines "
                    f"{', '.join(str(snp.line) for snp in found)}), where one is expected; only one "
                    f"of them is used",
                )
            )

    for node in sorted(set(by_node) - set(nodes)):
        findings.append(
            problem(
                MODIFY,
                "unknown-node",
                node,
                f"the SNP file annotates '{node}', which the tree file does not hold, so the SNP is not used",
            )
        )

    return findings


def check_alleles(scheme: Scheme) -> list[Finding]:
    """Whether each SNP can be called at all, and whether any two collide."""
    findings: list[Finding] = []
    positions: dict[tuple[str, str, int], list[str]] = {}

    for snp in scheme.checked_snps:
        alleles = {"derived": snp.derived, "ancestral": snp.ancestral}
        for label, base in alleles.items():
            if base not in BASES:
                findings.append(
                    problem(
                        MODIFY,
                        "allele",
                        snp.node,
                        snp.tag(
                            f"the {label} base is '{base}', not one of A, C, G and T. CanSNPer2 compares "
                            f"an aligned base against both columns, so a sample on this branch is placed "
                            f"above it without a warning"
                        ),
                    )
                )
        if snp.derived in BASES and snp.derived == snp.ancestral:
            findings.append(
                problem(
                    MODIFY,
                    "allele",
                    snp.node,
                    snp.tag(
                        f"the derived and ancestral base are both '{snp.derived}', so the branch cannot "
                        f"be told apart from its parent"
                    ),
                )
            )
        if snp.position < 1:
            findings.append(
                problem(
                    MODIFY,
                    "position",
                    snp.node,
                    snp.tag(f"the position is {snp.position}; positions are 1-based and positive"),
                )
            )

        positions.setdefault((snp.genome, snp.sequence, snp.position), []).append(snp.node)

    for (genome, sequence, position), nodes in sorted(positions.items()):
        if len(nodes) > 1:
            findings.append(
                problem(
                    MODIFY,
                    "shared-position",
                    ", ".join(sorted(nodes)),
                    f"{len(nodes)} nodes are annotated at {genome} {sequence}:{position}. One "
                    f"position cannot mark two branches, so at most one of them is right",
                )
            )

    return findings


def check_references(scheme: Scheme) -> list[Finding]:
    """The references file, and whether the SNPs name a genome it holds."""
    findings: list[Finding] = []
    known = scheme.references_by_genome

    for reference in scheme.references:
        node = reference.genome
        if "_" in reference.genome:
            findings.append(
                problem(
                    MODIFY,
                    "genome-name",
                    node,
                    f"the genome name '{reference.genome}' holds an underscore; CanSNPer2 recovers "
                    f"the reference name from an alignment file name by splitting on the first one",
                )
            )
        if not GENBANK.fullmatch(reference.genbank_id):
            findings.append(
                problem(
                    MODIFY,
                    "accession",
                    node,
                    f"'{reference.genbank_id}' is not a well formed GenBank accession, such as "
                    f"GCA_000007125.1; CanSNPer2-download builds NCBI paths from its digits",
                )
            )
        if not REFSEQ.fullmatch(reference.refseq_id):
            findings.append(
                problem(
                    MODIFY,
                    "accession",
                    node,
                    f"'{reference.refseq_id}' is not a well formed RefSeq accession, such as GCF_000007125.1",
                )
            )
        if not reference.assembly_name:
            findings.append(
                problem(
                    MODIFY,
                    "assembly-name",
                    node,
                    "the assembly name is empty, and the download needs it to name the file",
                )
            )

    for snp in scheme.checked_snps:
        if snp.genome not in known:
            findings.append(
                problem(
                    MODIFY,
                    "unknown-genome",
                    snp.node,
                    f"the SNP is annotated on genome '{snp.genome}', which the references file does not hold",
                )
            )

    used = {snp.genome for snp in scheme.snps}
    for genome in sorted(set(known) - used):
        findings.append(
            recommendation(
                MODIFY,
                "unused-genome",
                genome,
                "the references file holds this genome but no SNP is annotated on it",
            )
        )

    return findings


def resolve_sequences(scheme: Scheme, reference: Reference, aliases: Sequence[set[str]] | None) -> Scheme:
    """The scheme with each sequence named as the reference names it.

    A scheme exported from a reference with RefSeq names gives its sequences
    by their GenBank accessions, while the reference, the scores and the
    alignment made against it use the RefSeq names. A scheme sequence the
    reference does not hold is renamed to the one reference sequence its
    alias group, from --contig-aliases, holds. A group NCBI gives as not
    identical (Relationship <>) is not used, since a position on one
    sequence need not hold on the other. A name that cannot be resolved is
    kept, and the reference check reports it.
    """
    held = set(reference.names)
    renamed: dict[str, str] = {}
    for name in sorted({snp.sequence for snp in scheme.snps + scheme.snp_set} - held):
        found = {
            other
            for group in aliases or ()
            if name in group and not getattr(group, "differs", ())
            for other in group & held
        }
        if len(found) == 1:
            renamed[name] = found.pop()
    if not renamed:
        return scheme
    logging.info(
        "Read the scheme's sequences under the names of the reference: "
        + ", ".join(f"{name} as {contig}" for name, contig in renamed.items())
    )

    def named(snps):
        return [
            replace(snp, sequence=renamed[snp.sequence]) if snp.sequence in renamed else snp for snp in snps
        ]

    return replace(scheme, snps=named(scheme.snps), snp_set=named(scheme.snp_set))


def check_against_reference(scheme: Scheme, reference: Reference, genome: str | None = None) -> list[Finding]:
    """Whether the reference genome carries one of the two alleles at each position.

    This is the direct test for positions taken from a different version of the
    assembly: a scheme whose coordinates have moved reports a reference base
    that is neither the derived nor the ancestral allele at most of its SNPs.
    """
    findings: list[Finding] = []

    genomes = sorted({snp.genome for snp in scheme.snps})
    if genome is None:
        if len(genomes) == 1:
            genome = genomes[0]
        else:
            return [
                Finding(
                    RECOMMENDATION,
                    CHECK,
                    "reference-base",
                    "",
                    f"the scheme annotates {len(genomes)} reference genomes "
                    f"({', '.join(genomes)}), so --genome must say which one the reference file "
                    f"holds; the reference base was not checked",
                )
            ]

    checked = 0
    mismatched: list[str] = []
    derived_in_reference: list[str] = []
    unknown: dict[str, list[str]] = {}

    for snp in scheme.checked_snps:
        if snp.genome != genome or snp.position < 1:
            continue
        if snp.sequence not in reference.names:
            unknown.setdefault(snp.sequence, []).append(snp.node)
            continue
        contig = reference.contig(snp.sequence)
        if snp.position > contig.length:
            findings.append(
                problem(
                    MODIFY,
                    "position",
                    snp.node,
                    snp.tag(
                        f"the position {snp.position} is past the end of the reference sequence "
                        f"{contig.name}, which holds {contig.length} bases"
                    ),
                )
            )
            continue
        base = contig.base(snp.position - 1).upper()
        checked += 1
        if base == snp.derived.upper():
            derived_in_reference.append(snp.node)
        elif base != snp.ancestral.upper():
            mismatched.append(snp.node)
            findings.append(
                problem(
                    CHECK,
                    "reference-base",
                    snp.node,
                    snp.tag(
                        f"the reference carries '{base}' at {contig.name}:{snp.position}, which is "
                        f"neither the derived '{snp.derived}' nor the ancestral '{snp.ancestral}' "
                        f"base"
                    ),
                )
            )

    for sequence, nodes in sorted(unknown.items()):
        distinct = sorted(set(nodes))
        more = len(distinct) - 5
        held = len(reference) - 5
        findings.append(
            problem(
                CHECK,
                "sequence",
                ", ".join(distinct[:5]) + (f" and {more} more" if more > 0 else ""),
                f"{len(nodes)} SNP{'s are' if len(nodes) != 1 else ' is'} annotated on the "
                f"sequence '{sequence}', which the reference does not hold (it holds "
                f"{', '.join(reference.names[:5])}{f' and {held} more' if held > 0 else ''}), so "
                f"{'they were' if len(nodes) != 1 else 'it was'} not checked. Give the other names "
                f"of the reference's sequences with --contig-aliases <assembly_report.txt>, or use "
                f"a reference that names its sequences as the scheme does, as the GenBank assembly "
                f"of the scheme's references file does",
            )
        )

    if checked and len(mismatched) > checked / 2:
        findings.append(
            problem(
                MODIFY,
                "coordinates",
                "",
                f"{len(mismatched)} of {checked} SNPs disagree with the reference base. That is more "
                f"than half, so the positions were most likely taken against a different version of "
                f"the assembly and the whole scheme needs its coordinates remapped",
            )
        )

    if derived_in_reference:
        findings.append(
            Finding(
                RECOMMENDATION,
                CHECK,
                "reference-derived",
                "",
                f"{len(derived_in_reference)} SNPs at {len(set(derived_in_reference))} nodes have "
                f"the derived base in the reference genome itself, so the reference sits inside "
                f"those branches: {', '.join(sorted(set(derived_in_reference))[:5])}"
                f"{' and others' if len(set(derived_in_reference)) > 5 else ''}. "
                f"That is expected of a reference "
                f"chosen from within the scheme, and worth confirming otherwise",
            )
        )

    logging.info(f"Checked the reference base at {checked} of {len(scheme.checked_snps)} SNPs")
    return findings


def check_callability(scheme: Scheme, min_required_hits: int = MINIMUM_SNPS_PER_NODE) -> list[Finding]:
    """Whether CanSNPer2 can reach each node with its required number of hits.

    CanSNPer2 2.3.0 calls a node only when its path from the root carries at
    least ``--min_required_hits`` derived SNPs, three by default. A branch
    that lies closer to the root than that cannot be called however good its
    SNP is, which is what we met when typing the ported Brucella and Coxiella
    schemes. A node without a SNP is passed on the way to the nodes below it
    and does not count towards their SNPs, but cannot itself be the result.
    """
    findings: list[Finding] = []
    if scheme.hierarchy is None:
        return findings

    with_snp = set(scheme.snps_by_node)
    shallow: list[tuple[str, int]] = []
    callable_nodes = set()

    for node in scheme.hierarchy.nodes:
        if node not in with_snp:
            findings.append(
                problem(
                    EXTEND,
                    "unreachable",
                    node,
                    "the node has no SNP: CanSNPer2 passes it on the way to the nodes below it, "
                    "but it cannot be the result of typing, and it is not counted towards the "
                    "SNPs of their paths",
                )
            )
            continue
        snps = sum(1 for name in scheme.hierarchy.path(node) if name in with_snp)
        if snps < min_required_hits:
            shallow.append((node, snps))
        else:
            callable_nodes.add(node)

    if shallow:
        findings.append(
            recommendation(
                EXTEND,
                "min-required-hits",
                ", ".join(node for node, _ in shallow[:5]),
                f"{len(shallow)} nodes have fewer than {min_required_hits} SNPs on the path from the "
                f"root, so CanSNPer2 cannot call them with its default --min_required_hits of "
                f"{MINIMUM_SNPS_PER_NODE}. Either extend the scheme with canSNPs on the branches "
                f"above them, or type with a lower --min_required_hits and say so in the scheme's "
                f"documentation",
            )
        )

    if scheme.hierarchy.leaves and not callable_nodes & set(scheme.hierarchy.leaves):
        findings.append(
            problem(
                EXTEND,
                "min-required-hits",
                "",
                f"no leaf of the scheme has {min_required_hits} SNPs on its path from the root, so "
                f"with the default settings every sample is placed above the leaves",
            )
        )

    return findings


def check_spacing(scheme: Scheme, proximity: int = 50) -> list[Finding]:
    """SNPs of different nodes that lie close together.

    Only the SNPs of the SNP file are compared, those CanSNPer2 calls. The
    further ranks of a node's set are chosen apart from each other by export,
    and with several per node most positions have a neighbour of some node.
    """
    findings: list[Finding] = []

    by_sequence: dict[tuple[str, str], list[tuple[int, str]]] = {}
    for snp in scheme.snps:
        by_sequence.setdefault((snp.genome, snp.sequence), []).append((snp.position, snp.node))

    for (genome, sequence), entries in sorted(by_sequence.items()):
        entries.sort()
        for (first, left), (second, right) in itertools.pairwise(entries):
            gap = second - first
            if 0 < gap <= proximity:
                findings.append(
                    recommendation(
                        MODIFY,
                        "spacing",
                        f"{left}, {right}",
                        f"the two SNPs are {gap} bases apart on {genome} {sequence} ({first} and "
                        f"{second}). "
                        f"A deletion or a poorly aligned stretch there removes both from a genome "
                        f"at once, and each lies in the flank a PCR assay of the other would use",
                    )
                )

    return findings


def check_set_region(scheme: Scheme, region: int = 10000) -> list[Finding]:
    """Nodes most of whose ranked canSNPs lie in one stretch of the genome.

    The set of a node is meant to type a genome when some of its canSNPs
    lack a base. canSNPs close together lack it together: a deletion, an
    unassembled stretch or a poorly aligned region removes all of them from
    a genome at once, and the set then decides on the few left. A node is
    reported when more than half of its canSNPs, and at least three, lie
    within ``region`` bases of each other on one sequence.
    """
    findings: list[Finding] = []
    if region <= 0:
        return findings
    by_node: dict[tuple[str, str, str], list[int]] = {}
    for snp in scheme.checked_snps:
        by_node.setdefault((snp.node, snp.genome, snp.sequence), []).append(snp.position)
    counted: dict[str, int] = {}
    for (node, _genome, _sequence), positions in by_node.items():
        counted[node] = counted.get(node, 0) + len(positions)

    for (node, genome, sequence), positions in sorted(by_node.items()):
        total = counted[node]
        if total < 3:
            continue
        positions.sort()
        best, span = 0, (0, 0)
        start = 0
        for end, position in enumerate(positions):
            while position - positions[start] > region:
                start += 1
            if end - start + 1 > best:
                best, span = end - start + 1, (positions[start], position)
        if best >= 3 and best * 2 > total:
            findings.append(
                recommendation(
                    MODIFY,
                    "clustered-set",
                    node,
                    f"{best} of the node's {total} canSNPs lie within {span[1] - span[0]} bases on "
                    f"{genome} {sequence} ({span[0]}-{span[1]}). A deletion or a stretch without bases there "
                    f"removes them from a genome together, and the node is then decided by the "
                    f"rest; export with a larger --min-spacing spreads them further apart",
                )
            )
    return findings


def check_against_scores(
    scheme: Scheme, scores: Sequence[ScoredSnp], reference: Reference, genome: str | None = None
) -> list[Finding]:
    """Compare the scheme against scored candidates for the same nodes.

    Both the scores and the scheme give a position on one reference
    sequence, so they are matched by sequence, position and node. A
    candidate is offered only if it passed evaluate's filters, but the
    scheme's own SNP is looked up among all of them.
    The scores are against one reference genome, so only the scheme's SNPs on
    that genome are compared.
    """
    findings: list[Finding] = []
    annotated = sorted({snp.genome for snp in scheme.snps})
    if genome is None:
        if len(annotated) > 1:
            return [
                recommendation(
                    CHECK,
                    "not-scored",
                    "",
                    f"the scheme annotates {len(annotated)} reference genomes "
                    f"({', '.join(annotated)}), so --genome must say which one the scores are "
                    f"against; the scores were not compared",
                )
            ]
        genome = annotated[0] if annotated else None

    best: dict[str, ScoredSnp] = {}
    scored_at: dict[tuple[str, int, str], ScoredSnp] = {}
    for score in scores:
        if score.contig in reference.names:
            scored_at[(score.contig, score.position, score.node)] = score
    for score in scores:
        if not score.passed:
            continue
        if score.contig not in reference.names:
            continue
        if score.node not in best or score.total > best[score.node].total:
            best[score.node] = score

    for snp in scheme.snps:
        if snp.genome != genome:
            continue
        candidate = best.get(snp.node)
        current = scored_at.get((snp.sequence, snp.position, snp.node))

        if candidate is None:
            findings.append(
                recommendation(
                    CHECK,
                    "not-scored",
                    snp.node,
                    "the scores hold no canSNP for this node, so the scheme's SNP cannot be "
                    "compared and nothing better is on offer",
                )
            )
            continue
        if current is None:
            findings.append(
                recommendation(
                    CHECK,
                    "not-scored",
                    snp.node,
                    f"the scheme's position {snp.sequence}:{snp.position} is not among the scored "
                    f"candidates, so "
                    f"its quality is unknown; the best candidate for this node scores "
                    f"{candidate.total}",
                )
            )
            continue
        if candidate.total > current.total:
            findings.append(
                recommendation(
                    MODIFY,
                    "better-available",
                    snp.node,
                    f"the scheme's SNP at {snp.sequence}:{snp.position} scores {current.total}, "
                    f"while a candidate "
                    f"at {candidate.contig}:{candidate.position} scores {candidate.total}",
                )
            )

    scored_nodes = {score.node for score in scores if score.node and score.node != "None"}
    for node in sorted(scored_nodes - set(scheme.snps_by_node) - {""}):
        if node not in best:
            continue
        if scheme.hierarchy is not None and node not in scheme.hierarchy:
            continue
        findings.append(
            recommendation(
                EXTEND,
                "candidate-available",
                node,
                f"the node has no SNP in the scheme and a candidate scoring {best[node].total} is "
                f"available at {best[node].contig}:{best[node].position}",
            )
        )

    return findings


def check_masks(
    scheme: Scheme,
    reference: Reference,
    masks: Mapping[str, Sequence[tuple[int, int]]],
    genome: str | None = None,
) -> list[Finding]:
    """SNPs of the scheme that lie in a masked region, such as a recombinant one.

    A canSNP in a recombinant region may have been brought in by
    recombination rather than inherited down the branch, so it is the first
    to replace. The masks are given by reference sequence, as the scheme's
    positions are.
    """
    findings: list[Finding] = []
    for snp in scheme.checked_snps:
        if genome is not None and snp.genome != genome:
            continue
        if snp.sequence not in reference.names:
            continue  # reported by the reference check
        contig = reference.contig(snp.sequence)
        if not 1 <= snp.position <= contig.length:
            continue
        local = snp.position - 1
        for start, end in masks.get(contig.name, []):
            if start <= local <= end:
                findings.append(
                    recommendation(
                        MODIFY,
                        "masked",
                        snp.node,
                        snp.tag(
                            f"the SNP at {contig.name}:{snp.position} lies in a masked region "
                            f"({contig.name}:{start + 1}-{end + 1}), such as one recombination brought "
                            f"in; a SNP outside it would mark the branch more reliably"
                        ),
                    )
                )
                break
    return findings


def audit_scheme(
    scheme: Scheme,
    reference: Reference | None = None,
    scores: Sequence[ScoredSnp] | None = None,
    genome: str | None = None,
    min_required_hits: int = MINIMUM_SNPS_PER_NODE,
    proximity: int = 50,
    alignment: Mapping[str, str | ColumnRow] | None = None,
    clusters: Mapping[str, Sequence[str]] | None = None,
    min_clade_size: int = 2,
    min_clade_snps: int = 1,
    masks: Mapping[str, Sequence[tuple[int, int]]] | None = None,
    genome_conflicts: int = GENOME_CONFLICTS,
    region: int = 10000,
    new_clades: bool = True,
) -> list[Finding]:
    """Every check, in one list of findings, problems first.

    The checks against the genomes run when an alignment and the clusters
    placing its genomes are given, and need the reference genome as well.
    """
    findings: list[Finding] = [problem(CHECK, "unreadable", "", detail) for detail in scheme.problems]
    findings += check_structure(scheme)
    findings += check_alleles(scheme)
    findings += check_references(scheme)
    findings += check_callability(scheme, min_required_hits)
    findings += check_spacing(scheme, proximity)
    findings += check_set_region(scheme, region)

    if reference is not None:
        findings += check_against_reference(scheme, reference, genome)
        if scores is not None:
            findings += check_against_scores(scheme, scores, reference, genome)
        if masks:
            findings += check_masks(scheme, reference, masks, genome)
        if alignment is not None and clusters is not None:
            findings += check_against_alignment(
                scheme,
                alignment,
                reference,
                clusters,
                genome,
                min_clade_size,
                min_clade_snps,
                genome_conflicts,
                new_clades,
            )

    order = {PROBLEM: 0, RECOMMENDATION: 1}
    return sorted(findings, key=lambda finding: (order[finding.severity], finding.category, finding.node))


def _nodes_of(findings: Sequence[Finding], action: str) -> list[str]:
    """The nodes that findings of one action name, each once. A finding about
    two SNPs names both nodes in one field."""
    return sorted(
        {
            node.strip()
            for finding in findings
            if finding.action == action
            for node in finding.node.split(",")
            if node.strip()
        }
    )


def tally(findings: Sequence[Finding]) -> list[tuple[str, str, int]]:
    """How many findings there are of each category, problems first."""
    counts: dict[tuple[str, str], int] = {}
    for finding in findings:
        key = (finding.severity, finding.category)
        counts[key] = counts.get(key, 0) + 1
    order = {PROBLEM: 0, RECOMMENDATION: 1}
    return [
        (severity, category, found)
        for (severity, category), found in sorted(
            counts.items(), key=lambda item: (order[item[0][0]], -item[1], item[0][1])
        )
    ]


def summarize(scheme: Scheme, findings: Sequence[Finding]) -> list[str]:
    """The lines of the closing summary."""
    problems = [finding for finding in findings if finding.severity == PROBLEM]
    recommendations = [finding for finding in findings if finding.severity == RECOMMENDATION]
    to_extend = _nodes_of(findings, EXTEND)
    to_modify = _nodes_of(findings, MODIFY)

    further = sum(1 for snp in scheme.snp_set if snp.rank > 1)
    lines = [
        f"The '{scheme.organism}' scheme holds {len(scheme.nodes)} nodes, {len(scheme.snps)} "
        f"SNPs"
        + (f", {further} further ranked canSNPs" if further else "")
        + f" and {len(scheme.references)} reference genomes",
        f"{len(problems)} problems and {len(recommendations)} recommendations",
    ]
    for severity, label in ((PROBLEM, "Problems"), (RECOMMENDATION, "Recommendations")):
        counted = [f"{category} {found}" for kind, category, found in tally(findings) if kind == severity]
        if counted:
            lines.append(f"{label} by category: {', '.join(counted)}")
    if to_modify:
        lines.append(
            f"Nodes with a SNP to modify ({len(to_modify)}): "
            f"{', '.join(to_modify[:10])}" + (" and others" if len(to_modify) > 10 else "")
        )
    if to_extend:
        lines.append(
            f"Nodes where the scheme would be extended ({len(to_extend)}): "
            f"{', '.join(to_extend[:10])}" + (" and others" if len(to_extend) > 10 else "")
        )
    if not problems and not recommendations:
        lines.append("Nothing found that asks for a change")
    return lines


def write_findings(findings: Sequence[Finding], path: str) -> None:
    """Write the findings as a tab separated file."""
    with open(path, "w", newline="\n") as handle:
        handle.write("\t".join(Finding.HEADER) + "\n")
        handle.writelines("\t".join(finding.row()) + "\n" for finding in findings)
    logging.info(f"Wrote {len(findings)} findings to {path}")


def _scheme_columns(
    scheme: Scheme, reference_row: str, reference: Reference | None, genome: str | None
) -> set[int]:
    """The alignment columns of the scheme's canSNPs, from the alignment's first row.

    None are given when the scheme or the row cannot be mapped; the checks
    that follow report why.
    """
    if scheme.hierarchy is None or reference is None:
        return set()
    try:
        mapping = AlignmentMap(reference_row, reference)
    except ReferenceError:
        return set()
    return {
        column
        for snps in scheme_columns(scheme, mapping, genome).values()
        for _snp, column in snps
        if column is not None
    }


def _placement_problems(
    by_clade: bool,
    alignment: str | None,
    genome_map: str | None,
    tree: str | None,
    members: str | None,
    outgroups: Sequence[str],
    no_new_clades: bool,
    typing_report: str | None,
) -> list[str]:
    """The options that place the genomes of the alignment, and those that need one, refused in this order."""
    refusals = [
        (
            by_clade and not (tree and genome_map and outgroups),
            "--rule clade needs --tree, --genome-map with the scheme's genomes, and --outgroup",
        ),
        (outgroups and not by_clade, "--outgroup applies to --rule clade"),
        (by_clade and members, "--members applies to --rule jaccard"),
        (
            ((genome_map or tree) and not alignment) or (genome_map and tree and not by_clade),
            "checking the scheme against the genomes needs --alignment, and either "
            "--genome-map or --tree to say where the genomes sit in the scheme, or the "
            "scheme's own genomes file",
        ),
        (members and not tree, "--members applies to the leaves of --tree"),
        (no_new_clades and not alignment, "--no-new-clades applies to --alignment"),
        (
            typing_report and not alignment,
            "--typing-report types the genomes of --alignment, so it needs one",
        ),
    ]
    return [message for refused, message in refusals if refused]


def _option_problems(
    reference: str | None,
    scores: str | None,
    masks: str | None,
    masks_concatenated: bool,
    alignment: str | None,
    genome: str | None,
    placement_problems: Sequence[str],
) -> list[str]:
    """The option combinations the audit refuses, in the order they are checked.

    ``placement_problems`` are those of :func:`_placement_problems`, which come
    between the options that need the reference and the last two. Only the
    first problem is reported, so the order decides which one a user sees.
    """
    first = [
        (scores and not reference, "--scores needs --reference, the genome the scores were made against"),
        (masks and not reference, "--masks needs --reference, which names the sequences the masks lie on"),
        (alignment and not reference, "--alignment needs --reference, the genome the alignment is against"),
    ]
    last = [
        (
            genome and not reference,
            "--genome names the scheme's reference genome that --reference holds, so it needs --reference",
        ),
        (masks_concatenated and not masks, "--masks-concatenated applies to --masks"),
    ]
    return [
        *(message for refused, message in first if refused),
        *placement_problems,
        *(message for refused, message in last if refused),
    ]


@dataclass
class _Inputs:
    """The scheme and the files read beside it."""

    scheme: Scheme
    annotated: list[str]
    """The reference genomes the scheme's SNPs are annotated on."""
    sequence: Reference | None
    scored: list[ScoredSnp] | None
    masked: dict[str, list[tuple[int, int]]] | None


def _read_inputs(
    scheme_dir: str,
    organism: str | None,
    reference: str | None,
    scores: str | None,
    genome: str | None,
    masks: str | None,
    masks_concatenated: bool,
    contig_aliases: str | None,
) -> _Inputs:
    """Read the scheme, and the reference, scores and masks when given.

    Raises:
        OSError, ValueError, SchemeError, ReferenceError: If a file cannot be
            read, or the scheme annotates no SNP on ``genome``.
    """
    scheme = read_scheme(scheme_dir, organism)
    annotated = sorted({snp.genome for snp in scheme.snps})
    if genome and genome not in annotated:
        raise ValueError(
            f"the scheme annotates no SNP on genome '{genome}'; it names {', '.join(annotated) or 'none'}"
        )
    sequence = reference_parser(reference) if reference else None
    scored = vcf.read_scores(scores, sequence) if scores else None
    aliases = read_contig_aliases(contig_aliases) if contig_aliases else None
    masked = (
        read_masks(masks, sequence, masks_concatenated, aliases) if masks and sequence is not None else None
    )
    if sequence is not None:
        scheme = resolve_sequences(scheme, sequence, aliases)
    return _Inputs(scheme, annotated, sequence, scored, masked)


@dataclass
class _Placed:
    """The alignment as read, and where its genomes sit in the scheme; empty
    when no alignment is given."""

    aligned: Mapping[str, str | ColumnRow] | None = None
    clusters: dict[str, list[str]] | None = None
    findings: list[Finding] = field(default_factory=list)
    """The findings of a tree placement."""


def _read_aligned(
    alignment: str, scheme: Scheme, sequence: Reference, genome: str | None, every_column: bool
) -> tuple[Mapping[str, str | ColumnRow], dict[str, str] | None]:
    """The alignment, and the same alignment whole when every column was read."""
    if every_column:
        logging.info("Reading every column of the alignment")
        whole = read_alignment(alignment)
        return whole, whole
    logging.info("Reading only the scheme's columns of the alignment")
    aligned = read_alignment_columns(alignment, lambda row: _scheme_columns(scheme, row, sequence, genome))
    return aligned, None


def _clusters_from_genomes_file(
    scheme: Scheme, aligned: Mapping[str, str | ColumnRow], typing_report: str | None
) -> dict[str, list[str]] | None:
    """The genomes placed as the scheme's own genomes file lists them.

    The scheme lists the genomes it was made from; an outgroup sits at no
    node and counts as outside every node.

    Raises:
        ValueError: If the scheme lists no genomes and no typing report is
            asked for, so that there is nothing to do with the alignment.
    """
    if not scheme.genomes:
        if typing_report:
            logging.info(
                "No genome map, tree or genomes file says where the genomes sit, so "
                "they are typed with the scheme but its SNPs are not checked against "
                "a placement"
            )
            return None
        raise ValueError(
            f"the scheme lists no genomes ({scheme.organism}_genomes.txt), so "
            f"--genome-map or --tree must say where the genomes sit"
        )
    # Members are used where the alignment holds them, as one made with
    # RepGenR2's --all-genomes does, and passed over quietly where it does
    # not: they were never aligned.
    from cansnpscorer.helpers.names import resolve_names

    grouped: dict[str, list[str]] = {}
    listed = [entry for entry in scheme.genomes if entry[2] != "outgroup"]
    matched, _unmatched, _ambiguous = resolve_names(
        [accession for accession, _node, _role in listed], list(aligned)
    )
    skipped = 0
    for accession, node, role in listed:
        if role == "member" and accession not in matched:
            skipped += 1
            continue
        grouped.setdefault(node, []).append(accession)
    logging.info(
        f"Placed {sum(len(found) for found in grouped.values())} genomes as "
        f"the scheme's genomes file lists them"
        + (f"; {skipped} members are not in the alignment" if skipped else "")
    )
    return grouped


def _place_genomes(
    alignment: str,
    scheme: Scheme,
    sequence: Reference,
    genome: str | None,
    by_clade: bool,
    tree: str | None,
    genome_map: str | None,
    members: str | None,
    outgroups: Sequence[str],
    no_new_clades: bool,
    typing_report: str | None,
) -> _Placed:
    """Read the alignment and place its genomes in the scheme.

    Raises:
        OSError, ValueError, HierarchyError, ReferenceError, SchemeError: If
            a file cannot be read or the genomes cannot be placed.
    """
    placed = bool(by_clade or genome_map or tree or scheme.genomes)
    every_column = (placed and not no_new_clades) or bool(tree and not by_clade)
    aligned, whole = _read_aligned(alignment, scheme, sequence, genome, every_column)
    if by_clade:
        from cansnpscorer.cansnp_tree import place_clade

        # --rule clade needs both, which is checked before this is reached.
        assert tree and genome_map
        # The first record is the reference, which the tree need not hold.
        return _Placed(aligned, place_clade(scheme, tree, genome_map, outgroups, list(aligned)[1:]).clusters)
    if genome_map:
        return _Placed(aligned, read_genome_map(genome_map))
    if not tree:
        return _Placed(aligned, _clusters_from_genomes_file(scheme, aligned, typing_report))
    genome_tree = read_genome_tree(tree)
    folded = read_members(members, genome_tree) if members else None
    # A tree placement without the clade rule reads every column.
    assert whole is not None
    placement = place(scheme, genome_tree, whole, sequence, genome, members=folded)
    return _Placed(aligned, placement.clusters, placement.findings)


def _with_placement_findings(findings: list[Finding], placement_findings: Sequence[Finding]) -> list[Finding]:
    """The audit's findings with those of a tree placement, problems first."""
    if not placement_findings:
        return findings
    order = {PROBLEM: 0, RECOMMENDATION: 1}
    return sorted(
        [*findings, *placement_findings],
        key=lambda finding: (order[finding.severity], finding.category, finding.node),
    )


def _log_findings(scheme: Scheme, findings: Sequence[Finding]) -> None:
    """Log each finding, problems as errors, then the summary."""
    for finding in findings:
        if finding.severity == PROBLEM:
            logging.error(str(finding))
        else:
            logging.warning(str(finding))
    for line in summarize(scheme, findings):
        logging.info(line)


def _placed_by(
    scheme: Scheme, by_clade: bool, tree: str | None, genome_map: str | None, alignment: str | None
) -> str:
    """What placed the genomes, as the run table of the reports gives it."""
    if by_clade:
        return f"tree {tree}, by the clades of the genomes of {genome_map}"
    if genome_map:
        return f"genome map {genome_map}"
    if tree:
        return f"tree {tree}"
    if alignment and scheme.genomes:
        return "the scheme's genomes file"
    return NA


def _run_rows(
    scheme: Scheme,
    scheme_dir: str,
    reference: str | None,
    genomes: str,
    alignment: str | None,
    placed_by: str,
    scores: str | None,
    masks: str | None,
    min_required_hits: int,
    strictness: float,
    genome_conflicts: int,
) -> list[tuple[str, str]]:
    """The run table of the reports: what was read, and with which settings.

    ``genomes`` names the reference genomes of the scheme, as ``--genome``
    gives it or as the scheme annotates its SNPs.
    """
    return [
        ("cansnpcreator", __version__),
        ("scheme", f"{scheme.organism} in {scheme_dir}"),
        ("reference", str(reference) if reference else NA),
        ("reference genome of the scheme", genomes or NA),
        ("alignment", str(alignment) if alignment else NA),
        ("genomes placed by", placed_by),
        ("scores", str(scores) if scores else NA),
        ("masks", str(masks) if masks else NA),
        ("min_required_hits", str(min_required_hits)),
        ("strictness", f"{strictness:g}"),
        ("genome_conflicts", str(genome_conflicts)),
    ]


def _type_genomes(
    scheme: Scheme,
    placed: _Placed,
    sequence: Reference | None,
    genome: str | None,
    min_required_hits: int,
    genome_conflicts: int,
    strictness: float,
    typing_report: str | None,
    run: Sequence[tuple[str, str]],
) -> tuple[TypingReport | None, bool]:
    """Type the genomes of the alignment, and write the typing report when asked.

    Returns the typing, None when it could not be made, and whether a typing
    report that was asked for was not written.
    """
    typing = None
    try:
        # Typing follows an alignment that was read, and --alignment needs
        # --reference, which is checked before this is reached.
        assert placed.aligned is not None and sequence is not None
        typing = type_alignment(
            scheme,
            placed.aligned,
            sequence,
            genome,
            placed.clusters,
            min_required_hits,
            genome_conflicts,
            strictness=strictness,
            outgroups=[accession for accession, _node, role in scheme.genomes if role == "outgroup"],
        )
        if typing_report:
            write_report(typing, typing_report, run)
    except (OSError, TypingError, ReferenceError) as error:
        logging.error(f"No typing report was written: {error}")
        return typing, bool(typing_report)
    return typing, False


def _exit_status(findings: Sequence[Finding], strict: bool, not_typed: bool) -> int:
    """2 for a problem, or with ``strict`` any finding; 1 when the typing report failed; 0 otherwise."""
    if any(finding.severity == PROBLEM for finding in findings):
        return 2
    if strict and findings:
        return 2
    return 1 if not_typed else 0


def main(
    scheme_dir=None,
    organism=None,
    reference=None,
    scores=None,
    genome=None,
    min_required_hits=MINIMUM_SNPS_PER_NODE,
    proximity=50,
    output=None,
    strict=False,
    alignment=None,
    genome_map=None,
    min_clade_size=2,
    min_clade_snps=1,
    tree=None,
    members=None,
    masks=None,
    masks_concatenated=False,
    contig_aliases=None,
    genome_conflicts=GENOME_CONFLICTS,
    typing_report=None,
    report=None,
    region=10000,
    rule="jaccard",
    outgroups=(),
    no_new_clades=False,
    strictness=STRICTNESS,
) -> int:
    """Audit an existing scheme. Returns an exit status.

    With ``typing_report``, a prefix, every genome of the alignment is typed
    with the scheme and the tables of :mod:`cansnpscorer.typing_report` are
    written; that needs the alignment but no placement of its genomes. With
    ``report``, a path, the findings and the typing are also written as one
    Markdown file. With ``tree`` and ``rule`` 'clade', the genomes are
    placed by the clades of the scheme's genomes, which ``genome_map`` gives,
    in the tree rooted on ``outgroups``, as ``place --rule clade`` does.

    0 when nothing asks for a change, 1 when a file cannot be read, and 2 when
    the scheme holds a problem, or, with ``strict``, a recommendation.

    The checks for clades the scheme does not resolve read every column of
    the alignment. With ``no_new_clades`` they are skipped, and only the
    columns of the scheme's SNPs are read, which is what typing needs; that
    is also done when nothing places the genomes, so that only the typing
    runs. A tree placement reads every column, so it reads them all.
    """
    by_clade = rule == "clade"
    problems = _option_problems(
        reference,
        scores,
        masks,
        masks_concatenated,
        alignment,
        genome,
        _placement_problems(
            by_clade, alignment, genome_map, tree, members, outgroups, no_new_clades, typing_report
        ),
    )
    if problems:
        logging.error(problems[0])
        return 1
    try:
        inputs = _read_inputs(
            scheme_dir, organism, reference, scores, genome, masks, masks_concatenated, contig_aliases
        )
    except (OSError, ValueError, SchemeError, ReferenceError) as error:
        logging.error(error)
        return 1
    scheme, sequence = inputs.scheme, inputs.sequence

    placed = _Placed()
    if alignment:
        # --alignment needs --reference, which is checked above.
        assert sequence is not None
        try:
            placed = _place_genomes(
                alignment,
                scheme,
                sequence,
                genome,
                by_clade,
                tree,
                genome_map,
                members,
                outgroups,
                no_new_clades,
                typing_report,
            )
        except (OSError, ValueError, HierarchyError, ReferenceError, SchemeError) as error:
            logging.error(error)
            return 1

    findings = audit_scheme(
        scheme,
        sequence,
        inputs.scored,
        genome,
        min_required_hits,
        proximity,
        placed.aligned,
        placed.clusters,
        min_clade_size,
        min_clade_snps,
        inputs.masked,
        genome_conflicts,
        region,
        not no_new_clades,
    )
    findings = _with_placement_findings(findings, placed.findings)
    _log_findings(scheme, findings)

    if output:
        try:
            write_findings(findings, output)
        except OSError as error:
            logging.error(error)
            return 1

    # The findings are written first: a scheme that cannot be typed with, as
    # one with two roots, is exactly what they report.
    run = _run_rows(
        scheme,
        scheme_dir,
        reference,
        genome or ", ".join(inputs.annotated),
        alignment,
        _placed_by(scheme, by_clade, tree, genome_map, alignment),
        scores,
        masks,
        min_required_hits,
        strictness,
        genome_conflicts,
    )

    not_typed = False
    typing = None
    # --typing-report needs --alignment; a report is typed only when an
    # alignment was read.
    if placed.aligned is not None and (typing_report or report):
        typing, not_typed = _type_genomes(
            scheme,
            placed,
            sequence,
            genome,
            min_required_hits,
            genome_conflicts,
            strictness,
            typing_report,
            run,
        )

    if report:
        try:
            write_markdown(report, scheme, findings, typing, run)
            logging.info(f"Wrote the report to {report}")
        except OSError as error:
            logging.error(error)
            return 1

    return _exit_status(findings, strict, not_typed)
