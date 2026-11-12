"""Scoring of extracted canSNPs for their suitability as assay targets."""

import logging
import os
import shutil
from datetime import datetime
from typing import Dict, List, Optional, Sequence, Tuple

from cansnpscorer import vcf
from cansnpscorer.helpers.data_utils import apply_filters, summarize_scores
from cansnpscorer.helpers.file_utils import (
    read_gff3,
    read_hierarchy,
    reference_parser,
    trf_parser,
)
from cansnpscorer.helpers.scoring_methods import DEFAULT_WEIGHTS, contexts_for, score_snp
from cansnpscorer.hierarchy import HierarchyError
from cansnpscorer.model import Contig, Reference, ReferenceError
from cansnpscorer.validate import check_annotations, report

NUCLEOTIDES = frozenset("ACGT")


def annotations_by_contig(annotations: Dict[str, Sequence], reference: Reference,
                          description: str) -> Dict[str, Sequence]:
    """Attach annotations to the contigs of the reference, by exact sequence ID.

    GFF3, BED and VCF name a sequence by the ID the reference FASTA gives it,
    the first word of its header, so the names must agree. Annotations on a
    sequence the reference does not hold are left out with a warning.

    Raises:
        ValueError: If no annotated sequence is a contig of the reference,
            which means the file describes another genome or another naming.
    """
    if not annotations:
        logging.warning(f"The {description} file describes no sequences")
        return {}

    unmatched = sorted(set(annotations) - set(reference.names))
    if len(unmatched) == len(annotations):
        raise ValueError(
            f"none of the {len(annotations)} sequences of the {description} file "
            f"({', '.join(unmatched[:5])}) is a contig of the reference "
            f"({', '.join(reference.names[:5])}). Sequences are matched by ID, the first word "
            f"of each FASTA header, so the files must describe the same assembly under the "
            f"same names")
    if unmatched:
        logging.warning(
            f"The {description} file describes sequences that the reference does not hold, "
            f"and they were left out: {', '.join(unmatched[:5])}")
    missing = sorted(set(reference.names) - set(annotations))
    if missing:
        logging.warning(
            f"These contigs have no {description} annotation: {', '.join(missing[:5])}")
    return {name: annotations.get(name, []) for name in reference.names}


def locate_snps(snps: Sequence[vcf.Snp], reference: Reference,
                examples: int = 3) -> List[Tuple[Contig, int, str, Sequence[str]]]:
    """Place each SNP on its contig, as (contig, 0-based position, node, alleles).

    The alleles are (derived, ancestral) for a canSNP. A SNP that marks no node
    gets the node 'None'; it is kept, because it is still a neighbour of the
    canSNPs around it.

    Raises:
        ValueError: If a SNP falls outside its contig or disagrees with the
            reference base.
    """
    located: List[Tuple[Contig, int, str, Sequence[str]]] = []
    problems: List[str] = []
    for snp in snps:
        try:
            contig = reference.contig(snp.contig)
        except ReferenceError as error:
            raise ValueError(str(error)) from None
        local_position = snp.position - 1
        if not 0 <= local_position < contig.length:
            problems.append(f"{snp.contig}:{snp.position} is outside the contig")
            continue
        alleles = (snp.derived, snp.ancestral)
        base = contig.base(local_position)
        if base in NUCLEOTIDES and base not in alleles:
            problems.append(
                f"{contig.name}:{snp.position} is {base} in the reference but the SNP reports "
                f"{'/'.join(alleles)}")
            continue
        located.append((contig, local_position, snp.node or "None", alleles))
    if problems:
        raise ValueError(
            f"{len(problems)} of {len(snps)} SNPs do not fit the reference genome: "
            + "; ".join(problems[:examples]))
    return located


def back_up_existing(output: str) -> None:
    """Copy an output file out of the way before it is overwritten."""
    if not os.path.isfile(output):
        return
    timestamp = datetime.now().strftime('%Y%m%d%H%M%S')
    backup_name = f"{output}_{timestamp}.backup"
    shutil.copy(output, backup_name)
    logging.warning(f"{output} already existed; copied it to {backup_name}")


def main(snps=None, reference=None, genes=None, repeats=None, hierarchy=None, output=None,
         tandem_limit=50, gene_limit=0, gc_goal=0.5, cutoff=0, top_n=None,
         weights: Optional[Dict[str, float]] = None) -> int:
    """Score the canSNPs extract found. Returns a process exit status."""
    weights = dict(DEFAULT_WEIGHTS) if weights is None else weights
    try:
        genome = reference_parser(reference)
        parsed = vcf.read_snps(snps, genome)
        model = read_gff3(genes)
        gene_intervals = annotations_by_contig(model.genes, genome, "gene")
        cds = annotations_by_contig(model.cds, genome, "CDS") if model.cds else {}
        repeat_intervals = annotations_by_contig(trf_parser(repeats), genome, "repeat")
        tree = read_hierarchy(hierarchy) if hierarchy else None
    except (OSError, ValueError, HierarchyError, ReferenceError) as error:
        logging.error(error)
        return 1

    problems, warnings = check_annotations(gene_intervals, genome, "gene")
    found = check_annotations(repeat_intervals, genome, "repeat")
    if not report(problems + found[0], warnings + found[1], "annotations"):
        return 2

    try:
        located = locate_snps(parsed, genome)
    except (ValueError, ReferenceError) as error:
        logging.error(error)
        return 1

    if not cds:
        logging.warning(
            "The GFF3 file holds no CDS features, so no reading frame is known and every SNP "
            "scores 0 for the mutation term")
    # Neighbouring SNPs matter whether or not they are canonical, so the
    # contexts are built from every SNP and only then narrowed to the canonical
    # ones, which are the candidate assay targets.
    contexts = contexts_for(located, gene_intervals, repeat_intervals, cds)
    scores = [
        score_snp(context, tandem_limit, gene_limit, gc_goal, weights)
        for context in contexts
        if context.node != "None"
    ]
    apply_filters(scores, cutoff, top_n)
    passing = [score for score in scores if score.passed]
    logging.info(
        f"Scored {len(scores)} canSNPs out of {len(contexts)} SNPs; {len(passing)} pass the "
        f"cutoff of {cutoff}" + (f" and the best {top_n} per node" if top_n else ""))

    if tree is not None:
        summarize_scores(tree.nodes, passing)

    back_up_existing(output)
    written = vcf.write_scores(output, scores, genome, reference, cutoff, top_n, weights)
    logging.info(f"Wrote {written} scored canSNPs to {output}")
    return 0
