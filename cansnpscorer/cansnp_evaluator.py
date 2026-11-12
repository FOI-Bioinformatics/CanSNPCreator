"""Scoring of extracted canSNPs for their suitability as assay targets."""

import logging
import os
import shutil
from datetime import datetime
from typing import Dict, List, Optional, Sequence, Tuple

from cansnpscorer.helpers.data_utils import filter_scores_by_cutoff, summarize_scores
from cansnpscorer.helpers.file_utils import (
    ParsedSnp,
    read_gff3,
    read_hc_file,
    reference_parser,
    score_file_writer,
    snp_file_parser,
    trf_parser,
)
from cansnpscorer.helpers.scoring_methods import contexts_for, score_snp
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


def locate_snps(snps: Sequence[ParsedSnp], reference: Reference,
                examples: int = 3) -> List[Tuple[Contig, int, str, Sequence[str]]]:
    """Place each SNP on a contig and check it against the reference base.

    The reference base at a SNP must be one of the alleles reported for it.
    When it is not, the positions were produced against something other than
    this reference, and every score computed from them would be meaningless.

    Raises:
        ValueError: If any SNP falls outside the reference or disagrees with it.
    """
    located: List[Tuple[Contig, int, str, Sequence[str]]] = []
    outside: List[str] = []
    disagreeing: List[str] = []
    unknown = 0

    for snp in snps:
        node, alleles = snp.node, tuple(snp.alleles)

        if snp.contig is not None:
            try:
                contig = reference.contig(snp.contig)
            except ReferenceError as error:
                raise ValueError(str(error)) from None
            local_position = snp.position - 1
            if not 0 <= local_position < contig.length:
                outside.append(f"{snp.contig}:{snp.position}")
                continue
        else:
            if not 0 <= snp.global_position < reference.length:
                outside.append(str(snp.global_position + 1))
                continue
            contig, local_position = reference.locate(snp.global_position)

        base = contig.base(local_position)
        if base not in NUCLEOTIDES:
            unknown += 1
        elif base not in alleles:
            disagreeing.append(
                f"{contig.name}:{local_position + 1} is {base} in the reference but the "
                f"SNP reports {'/'.join(alleles)}")
        located.append((contig, local_position, node, alleles))

    if outside:
        raise ValueError(
            f"{len(outside)} SNPs fall outside the reference of {reference.length} bases "
            f"over {len(reference)} contigs, the first at {outside[0]}")
    if disagreeing:
        raise ValueError(
            f"{len(disagreeing)} of {len(snps)} SNPs disagree with the reference genome: "
            + "; ".join(disagreeing[:examples])
            + ". The SNP file was not produced against this reference")
    if unknown:
        logging.warning(f"{unknown} SNPs sit where the reference is not a plain base")

    return located


def back_up_existing(output: str) -> None:
    """Copy an output file out of the way before it is overwritten."""
    if not os.path.isfile(output):
        return
    timestamp = datetime.now().strftime('%Y%m%d%H%M%S')
    backup_name = f"{output}_{timestamp}.backup"
    shutil.copy(output, backup_name)
    logging.warning(f"{output} already existed; copied it to {backup_name}")


def main(snp_file=None, reference_file=None, gff_file=None, trf_file=None, hcfile=None,
         output_file=None, tandem_limit=50, gene_limit=0, gc_goal=0.5, cutoff_value=0,
         top_n_snp=None, weights: Optional[Dict[str, float]] = None) -> int:
    """Score the canSNPs in a SNP file. Returns a process exit status."""
    try:
        snps = snp_file_parser(snp_file)
        reference = reference_parser(reference_file)
        model = read_gff3(gff_file)
        genes = annotations_by_contig(model.genes, reference, "gene")
        cds = annotations_by_contig(model.cds, reference, "CDS") if model.cds else {}
        repeats = annotations_by_contig(trf_parser(trf_file), reference, "repeat")
        hierarchy = read_hc_file(hcfile)
    except (OSError, ValueError, HierarchyError, ReferenceError) as error:
        logging.error(error)
        return 1

    problems, warnings = check_annotations(genes, reference, "gene")
    found = check_annotations(repeats, reference, "repeat")
    if not report(problems + found[0], warnings + found[1], "annotations"):
        return 2

    try:
        located = locate_snps(snps, reference)
    except (ValueError, ReferenceError) as error:
        logging.error(error)
        return 1

    # Neighbouring SNPs matter whether or not they are canonical, so the
    # contexts are built from every SNP and only then narrowed to the canonical
    # ones, which are the candidate assay targets.
    if not cds:
        logging.warning(
            "The GFF3 file holds no CDS features, so no reading frame is known and every SNP "
            "scores 0 for the mutation term")
    contexts = contexts_for(located, genes, repeats, cds)
    scores = [
        score_snp(context, tandem_limit, gene_limit, gc_goal, weights)
        for context in contexts
        if context.node != "None"
    ]
    logging.info(f"Scored {len(scores)} canSNPs out of {len(contexts)} SNPs")

    filtered = filter_scores_by_cutoff(scores, cutoff_value, top_n_snp)
    logging.info(f"{len(filtered)} canSNPs passed the cutoff of {cutoff_value}")

    summarize_scores(hierarchy.nodes, filtered)

    back_up_existing(output_file)
    score_file_writer(filtered, output_file)
    return 0
