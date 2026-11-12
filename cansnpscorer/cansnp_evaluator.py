"""Scoring of extracted canSNPs for their suitability as assay targets."""

import logging
from bisect import bisect_right
from typing import Dict, List, Optional, Sequence, Set, Tuple

from cansnpscorer import vcf
from cansnpscorer.helpers.data_utils import apply_filters, summarize_scores
from cansnpscorer.helpers.file_utils import (
    merge_intervals,
    read_contig_aliases,
    read_gff3,
    read_hierarchy,
    read_intervals,
    reference_parser,
)
from cansnpscorer.helpers.scoring_methods import DEFAULT_WEIGHTS, contexts_for, score_snp
from cansnpscorer.hierarchy import HierarchyError
from cansnpscorer.model import Contig, Reference, ReferenceError
from cansnpscorer.validate import check_annotations, report

NUCLEOTIDES = frozenset("ACGT")


def contig_of(names: Sequence[str], reference: Reference,
              aliases: Optional[Sequence[Set[str]]] = None,
              sequences: Optional[Dict[str, str]] = None) -> Dict[str, str]:
    """Which contig of the reference each annotated sequence is.

    In order: the sequence's name is a contig ID; an alias of the name is,
    as NCBI's assembly report pairs a RefSeq accession with a GenBank one; or
    the annotation file carries the sequence itself, as Prokka and Bakta do
    after ##FASTA, and it is identical to a contig's. A name that matches
    nothing is left out of the result.
    """
    contigs = set(reference.names)
    by_sequence: Dict[str, List[str]] = {}
    if sequences:
        for contig in reference:
            by_sequence.setdefault(contig.sequence.upper(), []).append(contig.name)

    found: Dict[str, str] = {}
    for name in names:
        if name in contigs:
            found[name] = name
            continue
        for group in aliases or ():
            if name in group:
                matches = sorted(group & contigs)
                if len(matches) == 1:
                    found[name] = matches[0]
                    break
        if name in found:
            continue
        sequence = (sequences or {}).get(name)
        if sequence:
            matches = by_sequence.get(sequence.upper(), [])
            if len(matches) == 1:
                found[name] = matches[0]
    return found


def annotations_by_contig(annotations: Dict[str, Sequence], reference: Reference,
                          description: str,
                          aliases: Optional[Sequence[Set[str]]] = None,
                          sequences: Optional[Dict[str, str]] = None) -> Dict[str, Sequence]:
    """Attach annotations to the contigs of the reference.

    A sequence is matched by its ID, by an alias, or by its own sequence, as
    ``contig_of`` describes. Annotations on a sequence the reference does not
    hold are left out with a warning.

    Raises:
        ValueError: If no annotated sequence is a contig of the reference, or
            two annotated sequences are the same contig.
    """
    if not annotations:
        logging.warning(f"The {description} file describes no sequences")
        return {}

    found = contig_of(list(annotations), reference, aliases, sequences)
    unmatched = sorted(set(annotations) - set(found))
    if not found:
        raise ValueError(
            f"none of the {len(annotations)} sequences of the {description} file "
            f"({', '.join(unmatched[:5])}) is a contig of the reference "
            f"({', '.join(reference.names[:5])}). Sequences are matched by ID, the first word "
            f"of each FASTA header; a file from another naming of the same assembly needs "
            f"--contig-aliases, such as NCBI's assembly_report.txt, or the sequences after "
            f"##FASTA, as Prokka and Bakta write them")
    if unmatched:
        logging.warning(
            f"The {description} file describes sequences that the reference does not hold, "
            f"and they were left out: {', '.join(unmatched[:5])}")
    renamed = {name: contig for name, contig in found.items() if name != contig}
    if renamed:
        logging.info(
            f"Matched {len(renamed)} sequences of the {description} file to contigs by alias "
            f"or sequence: " + ", ".join(f"{a}={b}" for a, b in sorted(renamed.items())[:5]))

    attached: Dict[str, List] = {name: [] for name in reference.names}
    claimed: Dict[str, str] = {}
    for name, contig in found.items():
        if contig in claimed:
            raise ValueError(
                f"the {description} file describes the contig {contig} twice, as "
                f"{claimed[contig]} and as {name}")
        claimed[contig] = name
        attached[contig] = list(annotations[name])
    # A contig without annotation is reported by check_annotations where it
    # matters; a contig without repeats or masked regions is ordinary.
    return attached


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
        located.append((contig, local_position, snp.node or "None", alleles, snp.nearest))
    if problems:
        raise ValueError(
            f"{len(problems)} of {len(snps)} SNPs do not fit the reference genome: "
            + "; ".join(problems[:examples]))
    return located


def regions_over_contigs(regions: Dict[str, Sequence[Tuple[int, int]]],
                         reference: Reference) -> Dict[str, List[Tuple[int, int]]]:
    """Place regions given over the whole concatenated reference on its contigs.

    Gubbins, run on an alignment in reference coordinates such as RepGenR2's
    full_alignment.fasta, reports recombinant regions by alignment column,
    which is a position in the concatenated reference whatever sequence name
    the file gives. A region that spans two contigs is split between them.

    Raises:
        ValueError: If a region extends past the end of the reference.
    """
    placed: Dict[str, List[Tuple[int, int]]] = {name: [] for name in reference.names}
    for name, found in regions.items():
        for start, end in found:
            if end >= reference.length:
                raise ValueError(
                    f"the region {name}:{start + 1}-{end + 1} extends past the end of the "
                    f"reference, which holds {reference.length} bases over all its contigs")
            for contig in reference:
                if start <= contig.end - 1 and end >= contig.start:
                    placed[contig.name].append(
                        (max(start, contig.start) - contig.start, min(end, contig.end - 1) - contig.start))
    # Regions from several sequence names may overlap once placed; in_regions
    # needs them joined.
    return {name: merge_intervals(found) for name, found in placed.items()}


def read_masks(path: str, reference: Reference, concatenated: bool = False,
               aliases: Optional[Sequence[Set[str]]] = None
               ) -> Dict[str, List[Tuple[int, int]]]:
    """Masked regions by contig, from BED or GFF3 named by contig, or, with
    ``concatenated``, given over the whole concatenated reference."""
    regions = read_intervals(path)
    if concatenated:
        return regions_over_contigs(regions, reference)
    return annotations_by_contig(regions, reference, "mask", aliases)


def in_regions(position: int, regions: Sequence[Tuple[int, int]]) -> bool:
    """Whether a 0-based position lies in one of these sorted, joined regions."""
    index = bisect_right(regions, (position, float("inf"))) - 1
    return index >= 0 and regions[index][0] <= position <= regions[index][1]


def main(snps=None, reference=None, genes=None, repeats=None, hierarchy=None, output=None,
         tandem_limit=50, gene_limit=0, gc_goal=0.5, cutoff=0, top_n=None,
         weights: Optional[Dict[str, float]] = None, masks=None,
         masks_concatenated=False, contig_aliases=None, flank=50, gc_tolerance=0.1,
         spacing=1000) -> int:
    """Score the canSNPs extract found. Returns a process exit status."""
    weights = dict(DEFAULT_WEIGHTS) if weights is None else weights
    try:
        genome = reference_parser(reference)
        parsed = vcf.read_snps(snps, genome)
        aliases = read_contig_aliases(contig_aliases) if contig_aliases else None
        model = read_gff3(genes)
        gene_intervals = annotations_by_contig(model.genes, genome, "gene", aliases,
                                               model.sequences)
        cds = (annotations_by_contig(model.cds, genome, "CDS", aliases, model.sequences)
               if model.cds else {})
        repeat_intervals = annotations_by_contig(read_intervals(repeats), genome, "repeat",
                                                 aliases)
        masked = read_masks(masks, genome, masks_concatenated, aliases) if masks else {}
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
        return 2

    if not any(snp.nearest is not None for snp in parsed):
        logging.warning(
            "The SNPs carry no INFO/NEAREST, so the flanks are judged by the other SNPs of the "
            "file alone; extract records the variation among all the scheme's genomes")
    if not cds:
        logging.warning(
            "The GFF3 file holds no CDS features, so no reading frame is known and every SNP "
            "scores 0 for the mutation term")
    # Neighbouring SNPs matter whether or not they are canonical, so the
    # contexts are built from every SNP and only then narrowed to the canonical
    # ones, which are the candidate assay targets.
    contexts = contexts_for(located, gene_intervals, repeat_intervals, cds)
    scores = [
        score_snp(context, tandem_limit, gene_limit, gc_goal, weights, flank, gc_tolerance)
        for context in contexts
        if context.node != "None"
    ]
    # A canSNP in a masked region, such as a recombinant one, is marked before
    # the others are ranked, so that it does not take a place among the best.
    for score in scores:
        if in_regions(score.position - 1, masked.get(score.contig, [])):
            score.filters.append("Masked")
    if masked:
        logging.info(f"{sum('Masked' in s.filters for s in scores)} canSNPs lie in masked regions")
    apply_filters(scores, cutoff, top_n, spacing)
    passing = [score for score in scores if score.passed]
    logging.info(
        f"Scored {len(scores)} canSNPs out of {len(contexts)} SNPs; {len(passing)} pass the "
        f"cutoff of {cutoff}" + (f" and the best {top_n} per node, at least {spacing} bases "
                                 f"apart" if top_n else ""))

    if tree is not None:
        summarize_scores(tree.nodes, passing)

    written = vcf.write_scores(output, scores, genome, reference, cutoff, top_n, weights)
    logging.info(f"Wrote {written} scored canSNPs to {output}")
    return 0
