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

import logging
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from cansnpscorer import vcf
from cansnpscorer.alignment_audit import check_against_alignment
from cansnpscorer.cansnp_evaluator import read_masks
from cansnpscorer.cansnp_export import GENBANK, REFSEQ
from cansnpscorer.findings import (  # noqa: F401  re-exported for callers of this module
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
    read_contig_aliases,
    read_genome_map,
    reference_parser,
)
from cansnpscorer.helpers.scoring_methods import ScoredSnp
from cansnpscorer.hierarchy import HierarchyError
from cansnpscorer.model import Reference, ReferenceError
from cansnpscorer.scheme import Scheme, SchemeError, read_scheme
from cansnpscorer.tree_tools import place, read_genome_tree, read_members

BASES = frozenset("ACGT")



def check_structure(scheme: Scheme) -> List[Finding]:
    """The tree, and whether every node carries exactly one SNP."""
    findings: List[Finding] = []

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


def check_alleles(scheme: Scheme) -> List[Finding]:
    """Whether each SNP can be called at all, and whether any two collide."""
    findings: List[Finding] = []
    positions: Dict[Tuple[str, int], List[str]] = {}

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


def check_references(scheme: Scheme) -> List[Finding]:
    """The references file, and whether the SNPs name a genome it holds."""
    findings: List[Finding] = []
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
            f"{len(derived_in_reference)} SNPs have the derived base in the reference genome "
            f"itself, so the reference sits inside those branches: "
            f"{', '.join(sorted(derived_in_reference)[:5])}. That is expected of a reference "
            f"chosen from within the scheme, and worth confirming otherwise"))

    logging.info(f"Checked the reference base at {checked} of {len(scheme.checked_snps)} SNPs")
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
        findings.append(recommendation(
            EXTEND, "min-required-hits", ", ".join(node for node, _ in shallow[:5]),
            f"{len(shallow)} nodes have fewer than {min_required_hits} SNPs on the path from the "
            f"root, so CanSNPer2 cannot call them with its default --min_required_hits of "
            f"{MINIMUM_SNPS_PER_NODE}. Either extend the scheme with canSNPs on the branches "
            f"above them, or type with a lower --min_required_hits and say so in the scheme's "
            f"documentation"))

    if scheme.hierarchy.leaves and not callable_leaves:
        findings.append(problem(
            EXTEND, "min-required-hits", "",
            f"no leaf of the scheme has {min_required_hits} SNPs on its path from the root, so "
            f"with the default settings every sample is placed above the leaves"))

    return findings


def check_spacing(scheme: Scheme, proximity: int = 50) -> List[Finding]:
    """SNPs close enough together to interfere with one another's assay."""
    findings: List[Finding] = []

    by_genome: Dict[str, List[Tuple[int, str]]] = {}
    for snp in scheme.checked_snps:
        by_genome.setdefault(snp.genome, []).append((snp.position, snp.node))

    for genome, entries in sorted(by_genome.items()):
        entries.sort()
        for (first, left), (second, right) in zip(entries, entries[1:]):
            gap = second - first
            if 0 < gap <= proximity:
                findings.append(recommendation(
                    MODIFY, "spacing", f"{left}, {right}",
                    f"the two SNPs are {gap} bases apart on {genome} ({first} and {second}). A "
                    f"neighbouring SNP falls inside the primer or probe of an assay, so one of "
                    f"them is better replaced"))

    return findings


def check_against_scores(scheme: Scheme, scores: Sequence[ScoredSnp],
                         reference: Reference, genome: Optional[str] = None) -> List[Finding]:
    """Compare the scheme against scored candidates for the same nodes.

    The scores give positions by contig; the scheme counts them over the
    concatenated genome, so the candidates are converted before they are
    matched. A candidate is offered only if it passed evaluate's filters, but
    the scheme's own SNP is looked up among all of them, by position and node.
    The scores are against one reference genome, so only the scheme's SNPs on
    that genome are compared.
    """
    findings: List[Finding] = []
    annotated = sorted({snp.genome for snp in scheme.snps})
    if genome is None:
        if len(annotated) > 1:
            return [recommendation(
                CHECK, "not-scored", "",
                f"the scheme annotates {len(annotated)} reference genomes "
                f"({', '.join(annotated)}), so --genome must say which one the scores are "
                f"against; the scores were not compared")]
        genome = annotated[0] if annotated else None

    best: Dict[str, ScoredSnp] = {}
    scored_at: Dict[Tuple[int, str], ScoredSnp] = {}
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
                masks: Mapping[str, Sequence[Tuple[int, int]]],
                genome: Optional[str] = None) -> List[Finding]:
    """SNPs of the scheme that lie in a masked region, such as a recombinant one.

    A canSNP in a recombinant region may have been brought in by
    recombination rather than inherited down the branch, so it is the first
    to replace. Positions in the scheme count over the concatenated genome,
    so they are placed on their contig first.
    """
    findings: List[Finding] = []
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


def audit_scheme(scheme: Scheme, reference: Optional[Reference] = None,
                 scores: Optional[Sequence[ScoredSnp]] = None,
                 genome: Optional[str] = None,
                 min_required_hits: int = MINIMUM_SNPS_PER_NODE,
                 proximity: int = 50,
                 alignment: Optional[Mapping[str, str]] = None,
                 clusters: Optional[Mapping[str, Sequence[str]]] = None,
                 min_clade_size: int = 2, min_clade_snps: int = 1,
                 masks: Optional[Mapping[str, Sequence[Tuple[int, int]]]] = None
                 ) -> List[Finding]:
    """Every check, in one list of findings, problems first.

    The checks against the genomes run when an alignment and the clusters
    placing its genomes are given, and need the reference genome as well.
    """
    findings: List[Finding] = [
        problem(CHECK, "unreadable", "", detail) for detail in scheme.problems
    ]
    findings += check_structure(scheme)
    findings += check_alleles(scheme)
    findings += check_references(scheme)
    findings += check_callability(scheme, min_required_hits)
    findings += check_spacing(scheme, proximity)

    if reference is not None:
        findings += check_against_reference(scheme, reference, genome)
        if scores is not None:
            findings += check_against_scores(scheme, scores, reference, genome)
        if masks:
            findings += check_masks(scheme, reference, masks, genome)
        if alignment is not None and clusters is not None:
            findings += check_against_alignment(scheme, alignment, reference, clusters, genome,
                                                min_clade_size, min_clade_snps)

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


def main(scheme_dir=None, organism=None, reference=None, scores=None, genome=None,
         min_required_hits=MINIMUM_SNPS_PER_NODE, proximity=50, output=None,
         strict=False, alignment=None, genome_map=None, min_clade_size=2,
         min_clade_snps=1, tree=None, members=None, masks=None,
         masks_concatenated=False, contig_aliases=None) -> int:
    """Audit an existing scheme. Returns an exit status.

    0 when nothing asks for a change, 1 when a file cannot be read, and 2 when
    the scheme holds a problem, or, with ``strict``, a recommendation.
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
    if (genome_map or tree) and not alignment or (genome_map and tree):
        logging.error(
            "checking the scheme against the genomes needs --alignment, and either "
            "--genome-map or --tree to say where the genomes sit in the scheme, or the "
            "scheme's own genomes file")
        return 1
    if members and not tree:
        logging.error("--members applies to the leaves of --tree")
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
    placement_findings: List[Finding] = []
    if alignment:
        try:
            aligned = read_alignment(alignment)
            if genome_map:
                clusters = read_genome_map(genome_map)
            elif not tree:
                # The scheme lists the genomes it was made from; an outgroup
                # sits at no node and counts as outside every node.
                if not scheme.genomes:
                    logging.error(
                        f"the scheme lists no genomes ({scheme.organism}_genomes.txt), so "
                        f"--genome-map or --tree must say where the genomes sit")
                    return 1
                # Members are used where the alignment holds them, as one made
                # with RepGenR2's --all-genomes does, and passed over quietly
                # where it does not: they were never aligned.
                from cansnpscorer.helpers.names import resolve_names

                clusters = {}
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
                logging.info(
                    f"Placed {sum(len(found) for found in clusters.values())} genomes as the "
                    f"scheme's genomes file lists them"
                    + (f"; {skipped} members are not in the alignment" if skipped else ""))
            else:
                genome_tree = read_genome_tree(tree)
                folded = read_members(members, genome_tree) if members else None
                placement = place(scheme, genome_tree, aligned, sequence, genome,
                                  members=folded)
                clusters = placement.clusters
                placement_findings = placement.findings
        except (OSError, ValueError, HierarchyError, ReferenceError) as error:
            logging.error(error)
            return 1

    findings = audit_scheme(scheme, sequence, scored, genome, min_required_hits, proximity,
                            aligned, clusters, min_clade_size, min_clade_snps, masked)
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

    if any(finding.severity == PROBLEM for finding in findings):
        return 2
    if strict and findings:
        return 2
    return 0
