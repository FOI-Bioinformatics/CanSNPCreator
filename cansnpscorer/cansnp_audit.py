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
    scheme_columns,
    type_alignment,
    write_report,
)

BASES = frozenset("ACGT")



def check_structure(scheme: Scheme) -> list[Finding]:
    """The tree, and whether every node carries exactly one SNP."""
    findings: list[Finding] = []

    for node, parents in scheme.conflicts:
        findings.append(problem(
            MODIFY, "two-parents", node,
            f"the node is placed under more than one parent ({', '.join(parents)}). The "
            f"CanSNPer2 tree reader does not refuse this: it creates a second node of the same "
            f"name and attaches the SNP to one of them"))

    by_node = scheme.snps_by_node
    nodes = scheme.nodes

    for node in nodes:
        if ";" in node:
            findings.append(problem(
                MODIFY, "semicolon", node,
                "the node name holds a semicolon. The CanSNPer2 tree reader turns tabs into "
                "semicolons before splitting, so the name is read as two columns"))

    for node in nodes:
        found = by_node.get(node, [])
        if not found:
            findings.append(problem(
                EXTEND, "no-snp", node,
                "the node is in the tree but has no SNP, so it can never be called. Every node "
                "of a CanSNPer2 scheme carries exactly one"))
        elif len(found) > 1:
            findings.append(problem(
                MODIFY, "repeated-snp", node,
                f"{len(found)} SNPs annotate the node (lines "
                f"{', '.join(str(snp.line) for snp in found)}), where one is expected; only one "
                f"of them is used"))

    for node in sorted(set(by_node) - set(nodes)):
        findings.append(problem(
            MODIFY, "unknown-node", node,
            f"the SNP file annotates '{node}', which the tree file does not hold, so the SNP is "
            f"not used"))

    return findings


def check_alleles(scheme: Scheme) -> list[Finding]:
    """Whether each SNP can be called at all, and whether any two collide."""
    findings: list[Finding] = []
    positions: dict[tuple[str, int], list[str]] = {}

    for snp in scheme.checked_snps:
        alleles = {"derived": snp.derived, "ancestral": snp.ancestral}
        for label, base in alleles.items():
            if base not in BASES:
                findings.append(problem(
                    MODIFY, "allele", snp.node, snp.tag(
                    f"the {label} base is '{base}', not one of A, C, G and T. CanSNPer2 compares "
                    f"an aligned base against both columns, so a sample on this branch is placed "
                    f"above it without a warning")))
        if snp.derived in BASES and snp.derived == snp.ancestral:
            findings.append(problem(
                MODIFY, "allele", snp.node, snp.tag(
                f"the derived and ancestral base are both '{snp.derived}', so the branch cannot "
                f"be told apart from its parent")))
        if snp.position < 1:
            findings.append(problem(
                MODIFY, "position", snp.node, snp.tag(
                f"the position is {snp.position}; positions are 1-based and positive")))

        positions.setdefault((snp.genome, snp.position), []).append(snp.node)

    for (genome, position), nodes in sorted(positions.items()):
        if len(nodes) > 1:
            findings.append(problem(
                MODIFY, "shared-position", ", ".join(sorted(nodes)),
                f"{len(nodes)} nodes are annotated at {genome}:{position}. One position cannot "
                f"mark two branches, so at most one of them is right"))

    return findings


def check_references(scheme: Scheme) -> list[Finding]:
    """The references file, and whether the SNPs name a genome it holds."""
    findings: list[Finding] = []
    known = scheme.references_by_genome

    for reference in scheme.references:
        node = reference.genome
        if "_" in reference.genome:
            findings.append(problem(
                MODIFY, "genome-name", node,
                f"the genome name '{reference.genome}' holds an underscore; CanSNPer2 recovers "
                f"the reference name from an alignment file name by splitting on the first one"))
        if not GENBANK.fullmatch(reference.genbank_id):
            findings.append(problem(
                MODIFY, "accession", node,
                f"'{reference.genbank_id}' is not a well formed GenBank accession, such as "
                f"GCA_000007125.1; CanSNPer2-download builds NCBI paths from its digits"))
        if not REFSEQ.fullmatch(reference.refseq_id):
            findings.append(problem(
                MODIFY, "accession", node,
                f"'{reference.refseq_id}' is not a well formed RefSeq accession, such as "
                f"GCF_000007125.1"))
        if not reference.assembly_name:
            findings.append(problem(
                MODIFY, "assembly-name", node,
                "the assembly name is empty, and the download needs it to name the file"))

    for snp in scheme.checked_snps:
        if snp.genome not in known:
            findings.append(problem(
                MODIFY, "unknown-genome", snp.node,
                f"the SNP is annotated on genome '{snp.genome}', which the references file does "
                f"not hold"))

    used = {snp.genome for snp in scheme.snps}
    for genome in sorted(set(known) - used):
        findings.append(recommendation(
            MODIFY, "unused-genome", genome,
            "the references file holds this genome but no SNP is annotated on it"))

    return findings


def check_against_reference(scheme: Scheme, reference: Reference,
                            genome: str | None = None) -> list[Finding]:
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
            return [Finding(
                RECOMMENDATION, CHECK, "reference-base", "",
                f"the scheme annotates {len(genomes)} reference genomes "
                f"({', '.join(genomes)}), so --genome must say which one the reference file "
                f"holds; the reference base was not checked")]

    checked = 0
    mismatched: list[str] = []
    derived_in_reference: list[str] = []

    for snp in scheme.checked_snps:
        if snp.genome != genome or snp.position < 1:
            continue
        if snp.position > reference.length:
            findings.append(problem(
                MODIFY, "position", snp.node, snp.tag(
                f"the position {snp.position} is past the end of the reference genome, which "
                f"holds {reference.length} bases")))
            continue
        contig, local = reference.locate(snp.position - 1)
        base = contig.base(local).upper()
        checked += 1
        if base == snp.derived.upper():
            derived_in_reference.append(snp.node)
        elif base != snp.ancestral.upper():
            mismatched.append(snp.node)
            findings.append(problem(
                CHECK, "reference-base", snp.node, snp.tag(
                f"the reference carries '{base}' at {contig.name}:{local + 1} (position "
                f"{snp.position} of the concatenated genome), which is neither the derived "
                f"'{snp.derived}' nor the ancestral '{snp.ancestral}' base")))

    if checked and len(mismatched) > checked / 2:
        findings.append(problem(
            MODIFY, "coordinates", "",
            f"{len(mismatched)} of {checked} SNPs disagree with the reference base. That is more "
            f"than half, so the positions were most likely taken against a different version of "
            f"the assembly and the whole scheme needs its coordinates remapped"))

    if derived_in_reference:
        findings.append(Finding(
            RECOMMENDATION, CHECK, "reference-derived", "",
            f"{len(derived_in_reference)} SNPs at {len(set(derived_in_reference))} nodes have "
            f"the derived base in the reference genome itself, so the reference sits inside "
            f"those branches: {', '.join(sorted(set(derived_in_reference))[:5])}"
            f"{' and others' if len(set(derived_in_reference)) > 5 else ''}. "
            f"That is expected of a reference "
            f"chosen from within the scheme, and worth confirming otherwise"))

    logging.info(f"Checked the reference base at {checked} of {len(scheme.checked_snps)} SNPs")
    return findings


def check_callability(scheme: Scheme, min_required_hits: int = MINIMUM_SNPS_PER_NODE
                      ) -> list[Finding]:
    """Whether CanSNPer2 can reach each node with its required number of hits.

    CanSNPer2 2.2.1 and 2.2.2 call a node only when its path from the root carries at
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
            findings.append(problem(
                EXTEND, "unreachable", node,
                "the node has no SNP: CanSNPer2 passes it on the way to the nodes below it, "
                "but it cannot be the result of typing, and it is not counted towards the "
                "SNPs of their paths"))
            continue
        snps = sum(1 for name in scheme.hierarchy.path(node) if name in with_snp)
        if snps < min_required_hits:
            shallow.append((node, snps))
        else:
            callable_nodes.add(node)

    if shallow:
        findings.append(recommendation(
            EXTEND, "min-required-hits", ", ".join(node for node, _ in shallow[:5]),
            f"{len(shallow)} nodes have fewer than {min_required_hits} SNPs on the path from the "
            f"root, so CanSNPer2 cannot call them with its default --min_required_hits of "
            f"{MINIMUM_SNPS_PER_NODE}. Either extend the scheme with canSNPs on the branches "
            f"above them, or type with a lower --min_required_hits and say so in the scheme's "
            f"documentation"))

    if scheme.hierarchy.leaves and not callable_nodes & set(scheme.hierarchy.leaves):
        findings.append(problem(
            EXTEND, "min-required-hits", "",
            f"no leaf of the scheme has {min_required_hits} SNPs on its path from the root, so "
            f"with the default settings every sample is placed above the leaves"))

    return findings


def check_spacing(scheme: Scheme, proximity: int = 50) -> list[Finding]:
    """SNPs of different nodes that lie close together.

    Only the SNPs of the SNP file are compared, those CanSNPer2 calls. The
    further ranks of a node's set are chosen apart from each other by export,
    and with several per node most positions have a neighbour of some node.
    """
    findings: list[Finding] = []

    by_genome: dict[str, list[tuple[int, str]]] = {}
    for snp in scheme.snps:
        by_genome.setdefault(snp.genome, []).append((snp.position, snp.node))

    for genome, entries in sorted(by_genome.items()):
        entries.sort()
        for (first, left), (second, right) in itertools.pairwise(entries):
            gap = second - first
            if 0 < gap <= proximity:
                findings.append(recommendation(
                    MODIFY, "spacing", f"{left}, {right}",
                    f"the two SNPs are {gap} bases apart on {genome} ({first} and {second}). "
                    f"A deletion or a poorly aligned stretch there removes both from a genome "
                    f"at once, and each lies in the flank a PCR assay of the other would use"))

    return findings


def check_set_region(scheme: Scheme, region: int = 10000) -> list[Finding]:
    """Nodes most of whose ranked canSNPs lie in one stretch of the genome.

    The set of a node is meant to type a genome when some of its canSNPs
    lack a base. canSNPs close together lack it together: a deletion, an
    unassembled stretch or a poorly aligned region removes all of them from
    a genome at once, and the set then decides on the few left. A node is
    reported when more than half of its canSNPs, and at least three, lie
    within ``region`` bases of each other on one genome.
    """
    findings: list[Finding] = []
    if region <= 0:
        return findings
    by_node: dict[tuple[str, str], list[int]] = {}
    for snp in scheme.checked_snps:
        by_node.setdefault((snp.node, snp.genome), []).append(snp.position)
    counted: dict[str, int] = {}
    for (node, _genome), positions in by_node.items():
        counted[node] = counted.get(node, 0) + len(positions)

    for (node, genome), positions in sorted(by_node.items()):
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
            findings.append(recommendation(
                MODIFY, "clustered-set", node,
                f"{best} of the node's {total} canSNPs lie within {span[1] - span[0]} bases on "
                f"{genome} ({span[0]}-{span[1]}). A deletion or a stretch without bases there "
                f"removes them from a genome together, and the node is then decided by the "
                f"rest; export with a larger --min-spacing spreads them further apart"))
    return findings


def check_against_scores(scheme: Scheme, scores: Sequence[ScoredSnp],
                         reference: Reference, genome: str | None = None) -> list[Finding]:
    """Compare the scheme against scored candidates for the same nodes.

    The scores give positions by contig; the scheme counts them over the
    concatenated genome, so the candidates are converted before they are
    matched. A candidate is offered only if it passed evaluate's filters, but
    the scheme's own SNP is looked up among all of them, by position and node.
    The scores are against one reference genome, so only the scheme's SNPs on
    that genome are compared.
    """
    findings: list[Finding] = []
    annotated = sorted({snp.genome for snp in scheme.snps})
    if genome is None:
        if len(annotated) > 1:
            return [recommendation(
                CHECK, "not-scored", "",
                f"the scheme annotates {len(annotated)} reference genomes "
                f"({', '.join(annotated)}), so --genome must say which one the scores are "
                f"against; the scores were not compared")]
        genome = annotated[0] if annotated else None

    best: dict[str, ScoredSnp] = {}
    scored_at: dict[tuple[int, str], ScoredSnp] = {}
    for score in scores:
        if score.contig in reference.names:
            position = reference.contig(score.contig).start + score.position
            scored_at[(position, score.node)] = score
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
        current = scored_at.get((snp.position, snp.node))

        if candidate is None:
            findings.append(recommendation(
                CHECK, "not-scored", snp.node,
                "the scores hold no canSNP for this node, so the scheme's SNP cannot be "
                "compared and nothing better is on offer"))
            continue
        if current is None:
            findings.append(recommendation(
                CHECK, "not-scored", snp.node,
                f"the scheme's position {snp.position} is not among the scored candidates, so "
                f"its quality is unknown; the best candidate for this node scores "
                f"{candidate.total}"))
            continue
        if candidate.total > current.total:
            findings.append(recommendation(
                MODIFY, "better-available", snp.node,
                f"the scheme's SNP at {snp.position} scores {current.total}, while a candidate "
                f"at {candidate.contig}:{candidate.position} scores {candidate.total}"))

    scored_nodes = {score.node for score in scores if score.node and score.node != "None"}
    for node in sorted(scored_nodes - set(scheme.snps_by_node) - {""}):
        if node not in best:
            continue
        if scheme.hierarchy is not None and node not in scheme.hierarchy:
            continue
        findings.append(recommendation(
            EXTEND, "candidate-available", node,
            f"the node has no SNP in the scheme and a candidate scoring {best[node].total} is "
            f"available at {best[node].contig}:{best[node].position}"))

    return findings


def check_masks(scheme: Scheme, reference: Reference,
                masks: Mapping[str, Sequence[tuple[int, int]]],
                genome: str | None = None) -> list[Finding]:
    """SNPs of the scheme that lie in a masked region, such as a recombinant one.

    A canSNP in a recombinant region may have been brought in by
    recombination rather than inherited down the branch, so it is the first
    to replace. Positions in the scheme count over the concatenated genome,
    so they are placed on their contig first.
    """
    findings: list[Finding] = []
    for snp in scheme.checked_snps:
        if genome is not None and snp.genome != genome:
            continue
        if not 1 <= snp.position <= reference.length:
            continue
        contig, local = reference.locate(snp.position - 1)
        for start, end in masks.get(contig.name, []):
            if start <= local <= end:
                findings.append(recommendation(
                    MODIFY, "masked", snp.node, snp.tag(
                    f"the SNP at {contig.name}:{local + 1} lies in a masked region "
                    f"({contig.name}:{start + 1}-{end + 1}), such as one recombination brought "
                    f"in; a SNP outside it would mark the branch more reliably")))
                break
    return findings


def audit_scheme(scheme: Scheme, reference: Reference | None = None,
                 scores: Sequence[ScoredSnp] | None = None,
                 genome: str | None = None,
                 min_required_hits: int = MINIMUM_SNPS_PER_NODE,
                 proximity: int = 50,
                 alignment: Mapping[str, str] | None = None,
                 clusters: Mapping[str, Sequence[str]] | None = None,
                 min_clade_size: int = 2, min_clade_snps: int = 1,
                 masks: Mapping[str, Sequence[tuple[int, int]]] | None = None,
                 genome_conflicts: int = GENOME_CONFLICTS,
                 region: int = 10000, new_clades: bool = True) -> list[Finding]:
    """Every check, in one list of findings, problems first.

    The checks against the genomes run when an alignment and the clusters
    placing its genomes are given, and need the reference genome as well.
    """
    findings: list[Finding] = [
        problem(CHECK, "unreadable", "", detail) for detail in scheme.problems
    ]
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
            findings += check_against_alignment(scheme, alignment, reference, clusters, genome,
                                                min_clade_size, min_clade_snps,
                                                genome_conflicts, new_clades)

    order = {PROBLEM: 0, RECOMMENDATION: 1}
    return sorted(findings, key=lambda finding: (order[finding.severity], finding.category,
                                                 finding.node))


def _nodes_of(findings: Sequence[Finding], action: str) -> list[str]:
    """The nodes that findings of one action name, each once. A finding about
    two SNPs names both nodes in one field."""
    return sorted({node.strip() for finding in findings if finding.action == action
                   for node in finding.node.split(",") if node.strip()})


def tally(findings: Sequence[Finding]) -> list[tuple[str, str, int]]:
    """How many findings there are of each category, problems first."""
    counts: dict[tuple[str, str], int] = {}
    for finding in findings:
        key = (finding.severity, finding.category)
        counts[key] = counts.get(key, 0) + 1
    order = {PROBLEM: 0, RECOMMENDATION: 1}
    return [(severity, category, found) for (severity, category), found in sorted(
        counts.items(), key=lambda item: (order[item[0][0]], -item[1], item[0][1]))]


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
        counted = [f"{category} {found}" for kind, category, found in tally(findings)
                   if kind == severity]
        if counted:
            lines.append(f"{label} by category: {', '.join(counted)}")
    if to_modify:
        lines.append(f"Nodes with a SNP to modify ({len(to_modify)}): "
                     f"{', '.join(to_modify[:10])}"
                     + (" and others" if len(to_modify) > 10 else ""))
    if to_extend:
        lines.append(f"Nodes where the scheme would be extended ({len(to_extend)}): "
                     f"{', '.join(to_extend[:10])}"
                     + (" and others" if len(to_extend) > 10 else ""))
    if not problems and not recommendations:
        lines.append("Nothing found that asks for a change")
    return lines


def write_findings(findings: Sequence[Finding], path: str) -> None:
    """Write the findings as a tab separated file."""
    with open(path, "w", newline="\n") as handle:
        handle.write("\t".join(Finding.HEADER) + "\n")
        handle.writelines("\t".join(finding.row()) + "\n" for finding in findings)
    logging.info(f"Wrote {len(findings)} findings to {path}")


def _scheme_columns(scheme: Scheme, reference_row: str, reference: Reference | None,
                    genome: str | None) -> set[int]:
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
    return {column for snps in scheme_columns(scheme, mapping, genome).values()
            for _snp, column in snps if column is not None}


def main(scheme_dir=None, organism=None, reference=None, scores=None, genome=None,
         min_required_hits=MINIMUM_SNPS_PER_NODE, proximity=50, output=None,
         strict=False, alignment=None, genome_map=None, min_clade_size=2,
         min_clade_snps=1, tree=None, members=None, masks=None,
         masks_concatenated=False, contig_aliases=None,
         genome_conflicts=GENOME_CONFLICTS, typing_report=None, report=None,
         region=10000, rule="jaccard", outgroups=(), no_new_clades=False,
         strictness=STRICTNESS) -> int:
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
    if scores and not reference:
        logging.error(
            "--scores needs --reference, because the scores give positions by contig and the "
            "scheme counts them over the concatenated genome")
        return 1
    if masks and not reference:
        logging.error("--masks needs --reference, to place the scheme's positions on contigs")
        return 1
    if alignment and not reference:
        logging.error("--alignment needs --reference, the genome the alignment is against")
        return 1
    by_clade = rule == "clade"
    if by_clade and not (tree and genome_map and outgroups):
        logging.error("--rule clade needs --tree, --genome-map with the scheme's genomes, and "
                      "--outgroup")
        return 1
    if outgroups and not by_clade:
        logging.error("--outgroup applies to --rule clade")
        return 1
    if by_clade and members:
        logging.error("--members applies to --rule jaccard")
        return 1
    if ((genome_map or tree) and not alignment) or (genome_map and tree and not by_clade):
        logging.error(
            "checking the scheme against the genomes needs --alignment, and either "
            "--genome-map or --tree to say where the genomes sit in the scheme, or the "
            "scheme's own genomes file")
        return 1
    if members and not tree:
        logging.error("--members applies to the leaves of --tree")
        return 1
    if no_new_clades and not alignment:
        logging.error("--no-new-clades applies to --alignment")
        return 1
    if typing_report and not alignment:
        logging.error("--typing-report types the genomes of --alignment, so it needs one")
        return 1

    if genome and not reference:
        logging.error("--genome names the scheme's reference genome that --reference holds, so "
                      "it needs --reference")
        return 1
    if masks_concatenated and not masks:
        logging.error("--masks-concatenated applies to --masks")
        return 1
    try:
        scheme = read_scheme(scheme_dir, organism)
        annotated = sorted({snp.genome for snp in scheme.snps})
        if genome and genome not in annotated:
            logging.error(f"the scheme annotates no SNP on genome '{genome}'; it names "
                          f"{', '.join(annotated) or 'none'}")
            return 1
        sequence = reference_parser(reference) if reference else None
        scored = vcf.read_scores(scores, sequence) if scores else None
        aliases = read_contig_aliases(contig_aliases) if contig_aliases else None
        masked = read_masks(masks, sequence, masks_concatenated, aliases) if masks else None
    except (OSError, ValueError, SchemeError, ReferenceError) as error:
        logging.error(error)
        return 1

    aligned = clusters = None
    placement_findings: list[Finding] = []
    if alignment:
        try:
            placed = bool(by_clade or genome_map or tree or scheme.genomes)
            if (placed and not no_new_clades) or (tree and not by_clade):
                logging.info("Reading every column of the alignment")
                aligned = read_alignment(alignment)
            else:
                logging.info("Reading only the scheme's columns of the alignment")
                aligned = read_alignment_columns(
                    alignment, lambda row: _scheme_columns(scheme, row, sequence, genome))
            if by_clade:
                from cansnpscorer.cansnp_tree import place_clade

                # The first record is the reference, which the tree need not hold.
                clusters = place_clade(scheme, tree, genome_map, outgroups,
                                       list(aligned)[1:]).clusters
            elif genome_map:
                clusters = read_genome_map(genome_map)
            elif not tree:
                # The scheme lists the genomes it was made from; an outgroup
                # sits at no node and counts as outside every node.
                if not scheme.genomes and typing_report:
                    logging.info(
                        "No genome map, tree or genomes file says where the genomes sit, so "
                        "they are typed with the scheme but its SNPs are not checked against "
                        "a placement")
                elif not scheme.genomes:
                    logging.error(
                        f"the scheme lists no genomes ({scheme.organism}_genomes.txt), so "
                        f"--genome-map or --tree must say where the genomes sit")
                    return 1
                # Members are used where the alignment holds them, as one made
                # with RepGenR2's --all-genomes does, and passed over quietly
                # where it does not: they were never aligned.
                from cansnpscorer.helpers.names import resolve_names

                clusters = {} if scheme.genomes else None
                listed = [entry for entry in scheme.genomes if entry[2] != "outgroup"]
                matched, _unmatched, _ambiguous = resolve_names(
                    [accession for accession, _node, _role in listed], list(aligned))
                skipped = 0
                # The loop names its own variables: genome is the --genome option.
                for accession, node, role in listed:
                    if role == "member" and accession not in matched:
                        skipped += 1
                        continue
                    clusters.setdefault(node, []).append(accession)
                if clusters is not None:
                    logging.info(
                        f"Placed {sum(len(found) for found in clusters.values())} genomes as "
                        f"the scheme's genomes file lists them"
                        + (f"; {skipped} members are not in the alignment" if skipped else ""))
            else:
                genome_tree = read_genome_tree(tree)
                folded = read_members(members, genome_tree) if members else None
                placement = place(scheme, genome_tree, aligned, sequence, genome,
                                  members=folded)
                clusters = placement.clusters
                placement_findings = placement.findings
        except (OSError, ValueError, HierarchyError, ReferenceError, SchemeError) as error:
            logging.error(error)
            return 1

    findings = audit_scheme(scheme, sequence, scored, genome, min_required_hits, proximity,
                            aligned, clusters, min_clade_size, min_clade_snps, masked,
                            genome_conflicts, region, not no_new_clades)
    if placement_findings:
        order = {PROBLEM: 0, RECOMMENDATION: 1}
        findings = sorted(findings + placement_findings,
                          key=lambda finding: (order[finding.severity], finding.category,
                                               finding.node))

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

    # The findings are written first: a scheme that cannot be typed with, as
    # one with two roots, is exactly what they report.
    if by_clade:
        placed_by = f"tree {tree}, by the clades of the genomes of {genome_map}"
    elif genome_map:
        placed_by = f"genome map {genome_map}"
    elif tree:
        placed_by = f"tree {tree}"
    elif alignment and scheme.genomes:
        placed_by = "the scheme's genomes file"
    else:
        placed_by = NA
    run = [
        ("cansnpcreator", __version__),
        ("scheme", f"{scheme.organism} in {scheme_dir}"),
        ("reference", str(reference) if reference else NA),
        ("reference genome of the scheme", genome or ", ".join(annotated) or NA),
        ("alignment", str(alignment) if alignment else NA),
        ("genomes placed by", placed_by),
        ("scores", str(scores) if scores else NA),
        ("masks", str(masks) if masks else NA),
        ("min_required_hits", str(min_required_hits)),
        ("strictness", f"{strictness:g}"),
        ("genome_conflicts", str(genome_conflicts)),
    ]

    not_typed = False
    typing = None
    if typing_report or (report and aligned is not None):
        try:
            typing = type_alignment(
                scheme, aligned, sequence, genome, clusters, min_required_hits,
                genome_conflicts, strictness=strictness,
                outgroups=[accession for accession, _node, role in scheme.genomes
                           if role == "outgroup"])
            if typing_report:
                write_report(typing, typing_report, run)
        except (OSError, TypingError, ReferenceError) as error:
            logging.error(f"No typing report was written: {error}")
            not_typed = bool(typing_report)

    if report:
        try:
            write_markdown(report, scheme, findings, typing, run)
            logging.info(f"Wrote the report to {report}")
        except OSError as error:
            logging.error(error)
            return 1

    if any(finding.severity == PROBLEM for finding in findings):
        return 2
    if strict and findings:
        return 2
    return 1 if not_typed else 0
