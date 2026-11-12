"""Audit of an existing CanSNPer2-data scheme.

The validation script of the CanSNPer2-data repository answers whether a scheme
is well formed: whether the three files parse, whether every node carries a SNP
and whether the accessions look right. It does not answer the question asked of
a published scheme in practice, which is whether the scheme is any good: which
branches cannot be called at all, which SNPs can never be called derived, which
positions look as though they were taken from a different genome version, and
where a better SNP is available.

This module answers that question. It reads a scheme back in, optionally beside
the reference genome it is annotated against and a score file from the evaluate
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

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from cansnpscorer.cansnp_export import GENBANK, REFSEQ
from cansnpscorer.helpers.data_utils import MINIMUM_SNPS_PER_NODE
from cansnpscorer.helpers.file_utils import read_score_file, reference_parser
from cansnpscorer.helpers.scoring_methods import ScoredSnp
from cansnpscorer.model import Reference, ReferenceError
from cansnpscorer.scheme import Scheme, SchemeError, read_scheme

BASES = frozenset("ACGT")

PROBLEM = "problem"
RECOMMENDATION = "recommendation"

MODIFY = "modify"
EXTEND = "extend"
CHECK = "check"


@dataclass(frozen=True)
class Finding:
    """One thing the audit has to say about a scheme."""

    HEADER = ("severity", "action", "category", "node", "detail")

    severity: str
    action: str
    category: str
    node: str
    detail: str

    def row(self) -> List[str]:
        return [self.severity, self.action, self.category, self.node, self.detail]

    def __str__(self) -> str:
        where = f" [{self.node}]" if self.node else ""
        return f"{self.action}{where}: {self.detail}"


def _problem(action: str, category: str, node: str, detail: str) -> Finding:
    return Finding(PROBLEM, action, category, node, detail)


def _recommendation(action: str, category: str, node: str, detail: str) -> Finding:
    return Finding(RECOMMENDATION, action, category, node, detail)


def check_structure(scheme: Scheme) -> List[Finding]:
    """The tree, and whether every node carries exactly one SNP."""
    findings: List[Finding] = []

    for node, parents in scheme.conflicts:
        findings.append(_problem(
            MODIFY, "two-parents", node,
            f"the node is placed under more than one parent ({', '.join(parents)}). The "
            f"CanSNPer2 tree reader does not refuse this: it creates a second node of the same "
            f"name and attaches the SNP to one of them"))

    by_node = scheme.snps_by_node
    nodes = scheme.nodes

    for node in nodes:
        if ";" in node:
            findings.append(_problem(
                MODIFY, "semicolon", node,
                "the node name holds a semicolon. The CanSNPer2 tree reader turns tabs into "
                "semicolons before splitting, so the name is read as two columns"))

    for node in nodes:
        found = by_node.get(node, [])
        if not found:
            findings.append(_problem(
                EXTEND, "no-snp", node,
                "the node is in the tree but has no SNP, so it can never be called. Every node "
                "of a CanSNPer2 scheme carries exactly one"))
        elif len(found) > 1:
            findings.append(_problem(
                MODIFY, "repeated-snp", node,
                f"{len(found)} SNPs annotate the node (lines "
                f"{', '.join(str(snp.line) for snp in found)}), where one is expected; only one "
                f"of them is used"))

    for node in sorted(set(by_node) - set(nodes)):
        findings.append(_problem(
            MODIFY, "unknown-node", node,
            f"the SNP file annotates '{node}', which the tree file does not hold, so the SNP is "
            f"not used"))

    return findings


def check_alleles(scheme: Scheme) -> List[Finding]:
    """Whether each SNP can be called at all, and whether any two collide."""
    findings: List[Finding] = []
    positions: Dict[Tuple[str, int], List[str]] = {}

    for snp in scheme.snps:
        alleles = {"derived": snp.derived, "ancestral": snp.ancestral}
        for label, base in alleles.items():
            if base not in BASES:
                findings.append(_problem(
                    MODIFY, "allele", snp.node,
                    f"the {label} base is '{base}', not one of A, C, G and T. CanSNPer2 compares "
                    f"an aligned base against both columns, so a sample on this branch is placed "
                    f"above it without a warning"))
        if snp.derived in BASES and snp.derived == snp.ancestral:
            findings.append(_problem(
                MODIFY, "allele", snp.node,
                f"the derived and ancestral base are both '{snp.derived}', so the branch cannot "
                f"be told apart from its parent"))
        if snp.position < 1:
            findings.append(_problem(
                MODIFY, "position", snp.node,
                f"the position is {snp.position}; positions are 1-based and positive"))

        positions.setdefault((snp.genome, snp.position), []).append(snp.node)

    for (genome, position), nodes in sorted(positions.items()):
        if len(nodes) > 1:
            findings.append(_problem(
                MODIFY, "shared-position", ", ".join(sorted(nodes)),
                f"{len(nodes)} nodes are annotated at {genome}:{position}. One position cannot "
                f"mark two branches, so at most one of them is right"))

    return findings


def check_references(scheme: Scheme) -> List[Finding]:
    """The references file, and whether the SNPs name a genome it holds."""
    findings: List[Finding] = []
    known = scheme.references_by_genome

    for reference in scheme.references:
        node = reference.genome
        if "_" in reference.genome:
            findings.append(_problem(
                MODIFY, "genome-name", node,
                f"the genome name '{reference.genome}' holds an underscore; CanSNPer2 recovers "
                f"the reference name from an alignment file name by splitting on the first one"))
        if not GENBANK.match(reference.genbank_id):
            findings.append(_problem(
                MODIFY, "accession", node,
                f"'{reference.genbank_id}' is not a well formed GenBank accession, such as "
                f"GCA_000007125.1; CanSNPer2-download builds NCBI paths from its digits"))
        if not REFSEQ.match(reference.refseq_id):
            findings.append(_problem(
                MODIFY, "accession", node,
                f"'{reference.refseq_id}' is not a well formed RefSeq accession, such as "
                f"GCF_000007125.1"))
        if not reference.assembly_name:
            findings.append(_problem(
                MODIFY, "assembly-name", node,
                "the assembly name is empty, and the download needs it to name the file"))

    for snp in scheme.snps:
        if snp.genome not in known:
            findings.append(_problem(
                MODIFY, "unknown-genome", snp.node,
                f"the SNP is annotated on genome '{snp.genome}', which the references file does "
                f"not hold"))

    used = {snp.genome for snp in scheme.snps}
    for genome in sorted(set(known) - used):
        findings.append(_recommendation(
            MODIFY, "unused-genome", genome,
            "the references file holds this genome but no SNP is annotated on it"))

    return findings


def check_against_reference(scheme: Scheme, reference: Reference,
                            genome: Optional[str] = None) -> List[Finding]:
    """Whether the reference genome carries one of the two alleles at each position.

    This is the direct test for positions taken from a different version of the
    assembly: a scheme whose coordinates have moved reports a reference base
    that is neither the derived nor the ancestral allele at most of its SNPs.
    """
    findings: List[Finding] = []

    genomes = sorted({snp.genome for snp in scheme.snps})
    if genome is None:
        if len(genomes) == 1:
            genome = genomes[0]
        else:
            return [Finding(
                RECOMMENDATION, CHECK, "reference-base", "",
                f"the scheme annotates {len(genomes)} reference genomes "
                f"({', '.join(genomes)}), so --genome must say which one the reference file "
                f"holds; the reference base was not checked")]

    checked = 0
    mismatched: List[str] = []
    derived_in_reference: List[str] = []

    for snp in scheme.snps:
        if snp.genome != genome or snp.position < 1:
            continue
        if snp.position > reference.length:
            findings.append(_problem(
                MODIFY, "position", snp.node,
                f"the position {snp.position} is past the end of the reference genome, which "
                f"holds {reference.length} bases"))
            continue
        contig, local = reference.locate(snp.position - 1)
        base = contig.base(local).upper()
        checked += 1
        if base == snp.derived.upper():
            derived_in_reference.append(snp.node)
        elif base != snp.ancestral.upper():
            mismatched.append(snp.node)
            findings.append(_problem(
                CHECK, "reference-base", snp.node,
                f"the reference carries '{base}' at {contig.name}:{local + 1} (position "
                f"{snp.position} of the concatenated genome), which is neither the derived "
                f"'{snp.derived}' nor the ancestral '{snp.ancestral}' base"))

    if checked and len(mismatched) > checked / 2:
        findings.append(_problem(
            MODIFY, "coordinates", "",
            f"{len(mismatched)} of {checked} SNPs disagree with the reference base. That is more "
            f"than half, so the positions were most likely taken against a different version of "
            f"the assembly and the whole scheme needs its coordinates remapped"))

    if derived_in_reference:
        findings.append(Finding(
            RECOMMENDATION, CHECK, "reference-derived", "",
            f"{len(derived_in_reference)} SNPs have the derived base in the reference genome "
            f"itself, so the reference sits inside those branches: "
            f"{', '.join(sorted(derived_in_reference)[:5])}. That is expected of a reference "
            f"chosen from within the scheme, and worth confirming otherwise"))

    logging.info(f"Checked the reference base at {checked} of {len(scheme.snps)} SNPs")
    return findings


def check_callability(scheme: Scheme, min_required_hits: int = MINIMUM_SNPS_PER_NODE
                      ) -> List[Finding]:
    """Whether CanSNPer2 can reach each node with its required number of hits.

    CanSNPer2 calls the deepest node whose path from the root carries at least
    ``--min_required_hits`` SNPs, three by default. A branch that lies closer
    to the root than that cannot be called however good its SNP is, which is
    what we met when typing the ported Brucella and Coxiella schemes.
    """
    findings: List[Finding] = []
    if scheme.hierarchy is None:
        return findings

    with_snp = set(scheme.snps_by_node)
    shallow: List[Tuple[str, int]] = []

    for node in scheme.hierarchy.nodes:
        hits = sum(1 for name in scheme.hierarchy.path(node) if name in with_snp)
        if hits < min_required_hits:
            shallow.append((node, hits))

    callable_leaves = [node for node in scheme.hierarchy.leaves
                       if sum(1 for name in scheme.hierarchy.path(node) if name in with_snp)
                       >= min_required_hits]

    if shallow:
        findings.append(_recommendation(
            EXTEND, "min-required-hits", ", ".join(node for node, _ in shallow[:5]),
            f"{len(shallow)} nodes have fewer than {min_required_hits} SNPs on the path from the "
            f"root, so CanSNPer2 cannot call them with its default --min_required_hits of "
            f"{MINIMUM_SNPS_PER_NODE}. Either extend the scheme with canSNPs on the branches "
            f"above them, or type with a lower --min_required_hits and say so in the scheme's "
            f"documentation"))

    if scheme.hierarchy.leaves and not callable_leaves:
        findings.append(_problem(
            EXTEND, "min-required-hits", "",
            f"no leaf of the scheme has {min_required_hits} SNPs on its path from the root, so "
            f"with the default settings every sample is placed above the leaves"))

    return findings


def check_spacing(scheme: Scheme, proximity: int = 50) -> List[Finding]:
    """SNPs close enough together to interfere with one another's assay."""
    findings: List[Finding] = []

    by_genome: Dict[str, List[Tuple[int, str]]] = {}
    for snp in scheme.snps:
        by_genome.setdefault(snp.genome, []).append((snp.position, snp.node))

    for genome, entries in sorted(by_genome.items()):
        entries.sort()
        for (first, left), (second, right) in zip(entries, entries[1:]):
            gap = second - first
            if 0 < gap <= proximity:
                findings.append(_recommendation(
                    MODIFY, "spacing", f"{left}, {right}",
                    f"the two SNPs are {gap} bases apart on {genome} ({first} and {second}). A "
                    f"neighbouring SNP falls inside the primer or probe of an assay, so one of "
                    f"them is better replaced"))

    return findings


def check_against_scores(scheme: Scheme, scores: Sequence[ScoredSnp],
                         reference: Reference) -> List[Finding]:
    """Compare the scheme against scored candidates for the same nodes.

    The score file gives positions by contig; the scheme counts them over the
    concatenated genome, so the candidates are converted before they are
    matched.
    """
    findings: List[Finding] = []

    best: Dict[str, ScoredSnp] = {}
    global_position: Dict[int, ScoredSnp] = {}
    for score in scores:
        if score.contig not in reference.names:
            continue
        position = reference.contig(score.contig).start + score.position
        global_position[position] = score
        if score.node not in best or score.total > best[score.node].total:
            best[score.node] = score

    for snp in scheme.snps:
        candidate = best.get(snp.node)
        current = global_position.get(snp.position)

        if candidate is None:
            findings.append(_recommendation(
                CHECK, "not-scored", snp.node,
                "the score file holds no canSNP for this node, so the scheme's SNP cannot be "
                "compared and nothing better is on offer"))
            continue
        if current is None:
            findings.append(_recommendation(
                CHECK, "not-scored", snp.node,
                f"the scheme's position {snp.position} is not among the scored candidates, so "
                f"its quality is unknown; the best candidate for this node scores "
                f"{candidate.total}"))
            continue
        if candidate.total > current.total:
            findings.append(_recommendation(
                MODIFY, "better-available", snp.node,
                f"the scheme's SNP at {snp.position} scores {current.total}, while a candidate "
                f"at {candidate.contig}:{candidate.position} scores {candidate.total}"))

    scored_nodes = {score.node for score in scores if score.node and score.node != "None"}
    for node in sorted(scored_nodes - set(scheme.snps_by_node) - {""}):
        if node not in best:
            continue
        if scheme.hierarchy is not None and node not in scheme.hierarchy:
            continue
        findings.append(_recommendation(
            EXTEND, "candidate-available", node,
            f"the node has no SNP in the scheme and a candidate scoring {best[node].total} is "
            f"available at {best[node].contig}:{best[node].position}"))

    return findings


def audit_scheme(scheme: Scheme, reference: Optional[Reference] = None,
                 scores: Optional[Sequence[ScoredSnp]] = None,
                 genome: Optional[str] = None,
                 min_required_hits: int = MINIMUM_SNPS_PER_NODE,
                 proximity: int = 50) -> List[Finding]:
    """Every check, in one list of findings, problems first."""
    findings: List[Finding] = [
        _problem(CHECK, "unreadable", "", detail) for detail in scheme.problems
    ]
    findings += check_structure(scheme)
    findings += check_alleles(scheme)
    findings += check_references(scheme)
    findings += check_callability(scheme, min_required_hits)
    findings += check_spacing(scheme, proximity)

    if reference is not None:
        findings += check_against_reference(scheme, reference, genome)
        if scores is not None:
            findings += check_against_scores(scheme, scores, reference)

    order = {PROBLEM: 0, RECOMMENDATION: 1}
    return sorted(findings, key=lambda finding: (order[finding.severity], finding.category,
                                                 finding.node))


def summarize(scheme: Scheme, findings: Sequence[Finding]) -> List[str]:
    """The lines of the closing summary."""
    problems = [finding for finding in findings if finding.severity == PROBLEM]
    recommendations = [finding for finding in findings if finding.severity == RECOMMENDATION]
    to_extend = sorted({finding.node for finding in findings
                        if finding.action == EXTEND and finding.node})
    to_modify = sorted({finding.node for finding in findings
                        if finding.action == MODIFY and finding.node})

    lines = [
        f"The '{scheme.organism}' scheme holds {len(scheme.nodes)} nodes, {len(scheme.snps)} "
        f"SNPs and {len(scheme.references)} reference genomes",
        f"{len(problems)} problems and {len(recommendations)} recommendations",
    ]
    if to_modify:
        lines.append(f"SNPs to modify: {', '.join(to_modify[:10])}"
                     + (" and others" if len(to_modify) > 10 else ""))
    if to_extend:
        lines.append(f"Nodes where the scheme would be extended: {', '.join(to_extend[:10])}"
                     + (" and others" if len(to_extend) > 10 else ""))
    if not problems and not recommendations:
        lines.append("Nothing found that asks for a change")
    return lines


def write_findings(findings: Sequence[Finding], path: str) -> None:
    """Write the findings as a tab separated file."""
    with open(path, "w", newline="\n") as handle:
        handle.write("\t".join(Finding.HEADER) + "\n")
        for finding in findings:
            handle.write("\t".join(finding.row()) + "\n")
    logging.info(f"Wrote {len(findings)} findings to {path}")


def main(scheme_dir=None, organism=None, reference_file=None, score_file=None, genome=None,
         min_required_hits=MINIMUM_SNPS_PER_NODE, proximity=50, output=None,
         strict=False) -> int:
    """Audit an existing scheme. Returns an exit status.

    0 when nothing asks for a change, 1 when a file cannot be read, and 2 when
    the scheme holds a problem, or, with ``strict``, a recommendation.
    """
    try:
        scheme = read_scheme(scheme_dir, organism)
        reference = reference_parser(reference_file) if reference_file else None
        scores = read_score_file(score_file) if score_file else None
    except (OSError, ValueError, SchemeError, ReferenceError) as error:
        logging.error(error)
        return 1

    if score_file and reference is None:
        logging.error(
            "--score-file needs --reference-file, because the score file gives positions by "
            "contig and the scheme counts them over the concatenated genome")
        return 1

    findings = audit_scheme(scheme, reference, scores, genome, min_required_hits, proximity)

    for finding in findings:
        if finding.severity == PROBLEM:
            logging.error(str(finding))
        else:
            logging.warning(str(finding))

    for line in summarize(scheme, findings):
        logging.info(line)

    if output:
        try:
            write_findings(findings, output)
        except OSError as error:
            logging.error(error)
            return 1

    if any(finding.severity == PROBLEM for finding in findings):
        return 2
    if strict and findings:
        return 2
    return 0
