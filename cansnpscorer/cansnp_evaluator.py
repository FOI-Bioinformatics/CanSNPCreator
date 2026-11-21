"""Scoring of extracted canSNPs for their suitability as assay targets."""

import logging
from bisect import bisect_right
from collections.abc import Sequence

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
from cansnpscorer.helpers.kmers import DEFAULT_K, copy_numbers
from cansnpscorer.helpers.scoring_methods import (
    DEFAULT_PROFILE,
    PROFILES,
    contexts_for,
    score_snp,
)
from cansnpscorer.hierarchy import HierarchyError
from cansnpscorer.model import Contig, Reference, ReferenceError
from cansnpscorer.typing_report import judge_by_narrowest_group, read_observed
from cansnpscorer.validate import check_annotations, report

NUCLEOTIDES = frozenset("ACGT")


def contig_of(names: Sequence[str], reference: Reference,
              aliases: Sequence[set[str]] | None = None,
              sequences: dict[str, str] | None = None) -> dict[str, str]:
    """Which contig of the reference each annotated sequence is.

    In order: the sequence's name is a contig ID; an alias of the name is,
    as NCBI's assembly report pairs a RefSeq accession with a GenBank one; or
    the annotation file carries the sequence itself, as Prokka and Bakta do
    after ##FASTA, and it is identical to a contig's. A name that matches
    nothing is left out of the result.
    """
    contigs = set(reference.names)
    by_sequence: dict[str, list[str]] = {}
    if sequences:
        for contig in reference:
            by_sequence.setdefault(contig.sequence.upper(), []).append(contig.name)

    found: dict[str, str] = {}
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


def wrap_origin(intervals: Sequence, length: int, cds: bool = False) -> tuple[list, int]:
    """Split the features that run past the end of a circular contig.

    GFF3 gives a feature that spans the origin of a circular sequence an end
    beyond the sequence's length, and NCBI writes such features for circular
    chromosomes and plasmids. Each is split into the part up to the contig's
    end and the part from its start. A CDS keeps its reading frame: the part
    read first keeps the phase, and the other gets the phase that continues
    the codons across the origin.

    Only a feature that starts within the contig, runs past its end, and is
    shorter than the contig is split, so that its two parts do not overlap;
    any other stays as it is, for check_annotations to report.

    Returns:
        The intervals, sorted by start and end, and how many were split.
    """
    kept = []
    split = 0
    for interval in intervals:
        start, end = interval[0], interval[1]
        if not (0 <= start < length <= end and end - length < start):
            kept.append(interval)
            continue
        split += 1
        upper, lower = (start, length - 1), (0, end - length)
        rest = tuple(interval[2:])
        if cds:
            strand, phase = interval[2], interval[3]
            # The part read first keeps the phase; the bases it carries past
            # its last whole codon decide where the other part's codon begins.
            first, second = (upper, lower) if strand != "-" else (lower, upper)
            carried = (first[1] - first[0] + 1 - phase) % 3
            following = (3 - carried) % 3
            kept.append((*first, strand, phase, *interval[4:]))
            kept.append((*second, strand, following, *interval[4:]))
        else:
            kept.append((*upper, *rest))
            kept.append((*lower, *rest))
    kept.sort(key=lambda interval: (interval[0], interval[1]))
    return kept, split


def annotations_by_contig(annotations: dict[str, Sequence], reference: Reference,
                          description: str,
                          aliases: Sequence[set[str]] | None = None,
                          sequences: dict[str, str] | None = None) -> dict[str, Sequence]:
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

    attached: dict[str, list] = {name: [] for name in reference.names}
    claimed: dict[str, str] = {}
    for name, contig in found.items():
        if contig in claimed:
            raise ValueError(
                f"the {description} file describes the contig {contig} twice, as "
                f"{claimed[contig]} and as {name}")
        claimed[contig] = name
        attached[contig], split = wrap_origin(annotations[name], reference.contig(contig).length,
                                              cds=description == "CDS")
        if split:
            logging.warning(
                f"{split} {description} intervals on '{contig}' run past its end, as GFF3 "
                f"gives a feature across the origin of a circular sequence; each was split at "
                f"the origin into two")
    # A contig without annotation is reported by check_annotations where it
    # matters; a contig without repeats or masked regions is ordinary.
    return attached


def locate_snps(snps: Sequence[vcf.Snp], reference: Reference,
                examples: int = 3) -> list[tuple[Contig, int, str, Sequence[str]]]:
    """Place each SNP on its contig, as (contig, 0-based position, node, alleles).

    The alleles are (derived, ancestral) for a canSNP. A SNP that marks no node
    gets the node 'None'; it is kept, because it is still a neighbour of the
    canSNPs around it.

    Raises:
        ValueError: If a SNP falls outside its contig or disagrees with the
            reference base.
    """
    located: list[tuple[Contig, int, str, Sequence[str]]] = []
    problems: list[str] = []
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


def locus_at(position: int, loci: Sequence[tuple[int, int]], starts: Sequence[int],
             longest: int) -> tuple[int, int] | None:
    """The gene that holds a position: of overlapping genes, the first by start.

    ``loci`` are sorted by start, ``starts`` are their starts and ``longest``
    the length of the longest, which bounds how far back a gene can begin.
    """
    found = None
    index = bisect_right(starts, position) - 1
    while index >= 0 and position - loci[index][0] < longest:
        if loci[index][1] >= position:
            found = loci[index]
        index -= 1
    return found


def regions_over_contigs(regions: dict[str, Sequence[tuple[int, int]]],
                         reference: Reference) -> dict[str, list[tuple[int, int]]]:
    """Place regions given over the whole concatenated reference on its contigs.

    Gubbins, run on an alignment in reference coordinates such as RepGenR2's
    full_alignment.fasta, reports recombinant regions by alignment column,
    which is a position in the concatenated reference whatever sequence name
    the file gives. A region that spans two contigs is split between them.

    Raises:
        ValueError: If a region extends past the end of the reference.
    """
    placed: dict[str, list[tuple[int, int]]] = {name: [] for name in reference.names}
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
               aliases: Sequence[set[str]] | None = None
               ) -> dict[str, list[tuple[int, int]]]:
    """Masked regions by contig, from BED or GFF3 named by contig, or, with
    ``concatenated``, given over the whole concatenated reference."""
    regions = read_intervals(path)
    if concatenated:
        return regions_over_contigs(regions, reference)
    return annotations_by_contig(regions, reference, "mask", aliases)


def in_regions(position: int, regions: Sequence[tuple[int, int]]) -> bool:
    """Whether a 0-based position lies in one of these sorted, joined regions."""
    index = bisect_right(regions, (position, float("inf"))) - 1
    return index >= 0 and regions[index][0] <= position <= regions[index][1]


def main(snps=None, reference=None, genes=None, repeats=None, hierarchy=None, output=None,
         tandem_limit=50, gene_limit=0, gc_goal=0.5, cutoff=0, top_n=None,
         weights: dict[str, float] | None = None, masks=None,
         masks_concatenated=False, contig_aliases=None, flank=50, gc_tolerance=0.1,
         spacing=1000, copy_kmer=DEFAULT_K, observed=None) -> int:
    """Score the canSNPs extract found. Returns a process exit status.

    ``cutoff`` None applies no cutoff. ``observed`` is the _snps.tsv table of
    a typing report of a scheme holding these canSNPs: those that disagreed
    with their node, or had no base, in its genomes are ranked after the
    others of their node.

    With ``copy_kmer`` above 0, the copy number of each canSNP's context is
    counted with k-mers of that size, and canSNPs in repeated sequence are
    ranked after the others of their node.
    """
    weights = dict(PROFILES[DEFAULT_PROFILE][0]) if weights is None else weights
    try:
        genome = reference_parser(reference)
        vcf.check_contig_names(genome)
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
        behaviour = read_observed(observed) if observed else {}
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
    elif flank > vcf.NEAREST_CAP:
        # NEAREST is a lower bound at the cap, so the neighbour term is one too.
        logging.warning(
            f"--flank {flank} is larger than the {vcf.NEAREST_CAP} bases up to which extract "
            f"records the distance to the nearest variation (INFO/NEAREST), so variation "
            f"further away is not known: a canSNP with none within {vcf.NEAREST_CAP} bases "
            f"gets {vcf.NEAREST_CAP / flank:.2f} of the full neighbour term, and distances "
            f"beyond {vcf.NEAREST_CAP} are not told apart")
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
    missing = {(snp.contig, snp.position): snp.missing for snp in parsed}
    for score in scores:
        score.missing = missing.get((score.contig, score.position), 0)
    if copy_kmer:
        copies = copy_numbers(genome, [(score.contig, score.position - 1) for score in scores],
                              copy_kmer)
        for score in scores:
            score.copies = copies[(score.contig, score.position - 1)]
        repeated = sum(1 for score in scores if score.copies > 1)
        logging.info(f"{repeated} canSNPs lie in sequence that occurs more than once in the "
                     f"reference ({copy_kmer}-mers, both strands); within a node they are "
                     f"ranked after those in unique sequence")
    # The gene each canSNP lies in, from the genes one by one: joined intervals
    # would make a run of abutting genes one locus.
    contig_names = contig_of(list(model.loci), genome, aliases, model.sequences)
    loci = {contig_names[name]: wrap_origin(found, genome.contig(contig_names[name]).length)[0]
            for name, found in model.loci.items() if name in contig_names}
    lookup = {contig: ([start for start, _end in found],
                       max(end - start + 1 for start, end in found))
              for contig, found in loci.items() if found}
    for score in scores:
        if score.contig not in lookup:
            continue
        starts, longest = lookup[score.contig]
        gene = locus_at(score.position - 1, loci[score.contig], starts, longest)
        if gene is not None:
            score.locus = f"{gene[0] + 1}-{gene[1] + 1}"
    if behaviour:
        behaviour, regrouped = judge_by_narrowest_group(behaviour)
        if regrouped:
            logging.warning(
                f"At {len(regrouped)} nodes the canSNPs are derived in groups of genomes of "
                f"different width and no group holds more than half of them, so the genomes "
                f"typed join these branches partway. Their canSNPs were judged against the "
                f"narrowest group and not against the majority of the node; representatives "
                f"of the other groups would divide these branches: {', '.join(regrouped)}")
        matched = 0
        other_alleles = 0
        # The report names each sequence as the scheme does, which for an
        # exported scheme is its GenBank accession; the scores name it as the
        # reference does, which the alias table may relate to it.
        named: dict[str, list[str]] = {}
        for score in scores:
            if score.contig not in named:
                # A pair NCBI gives as not identical (<>) is not related.
                named[score.contig] = [score.contig, *sorted(
                    {name for group in aliases or ()
                     if score.contig in group and not getattr(group, "differs", ())
                     for name in group} - {score.contig})]
        for score in scores:
            seen = next((behaviour[key] for key in (
                (score.node, name, score.position) for name in named[score.contig])
                if key in behaviour), None)
            if seen is not None and (seen.derived, seen.ancestral) != (
                    (score.derived or "").upper(), (score.ancestral or "").upper()):
                # The same node and position with other alleles is another
                # SNP: a report made against another reference genome.
                other_alleles += 1
                seen = None
            if seen is not None:
                matched += 1
                score.discordant, score.observed_missing = seen.discordant, seen.missing
                score.observed_genomes = seen.genomes
        logging.info(
            f"{matched} of {len(scores)} canSNPs were observed in the typing report, "
            f"{sum(1 for s in scores if s.discordant)} of them disagreeing with their node in "
            f"some genome and {sum(1 for s in scores if s.observed_missing)} lacking a base in "
            f"some; within a node these are ranked after the others")
        if other_alleles:
            logging.warning(
                f"{other_alleles} rows of {observed} name the node and position of a canSNP "
                f"here with other alleles, and were not used; the report may have been made "
                f"against another reference genome")
        if 0 < matched < len(scores):
            logging.warning(
                f"{len(scores) - matched} canSNPs were not in the scheme the report typed "
                f"with. Within a node they are ranked after the observed ones that held and "
                f"before those that did not; for a ranking of every candidate, type with a "
                f"scheme exported with --snps-per-node all")
        if not matched:
            logging.warning(
                f"No canSNP of {observed} is among those scored; the table should come from "
                f"a scheme exported from these canSNPs against the same reference")
    incomplete = sum(1 for score in scores if score.missing)
    if incomplete:
        logging.info(f"{incomplete} canSNPs lack a base in some genomes; within a node they "
                     f"are ranked after those called in every genome")
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
        f"Scored {len(scores)} canSNPs out of {len(contexts)} SNPs; {len(passing)} pass "
        + (f"the cutoff of {cutoff}" if cutoff is not None else "with no cutoff") + (f" and the best {top_n} per node, at least {spacing} bases "
                                 f"apart" if top_n else ""))

    if tree is not None:
        summarize_scores(tree.nodes, passing)

    written = vcf.write_scores(output, scores, genome, reference, cutoff, top_n, weights,
                               copies=bool(copy_kmer), observed=bool(behaviour),
                               spacing=spacing)
    logging.info(f"Wrote {written} scored canSNPs to {output}")
    return 0
