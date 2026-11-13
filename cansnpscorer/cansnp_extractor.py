"""Extraction of canSNPs from a whole genome alignment."""

import logging
import os
from datetime import datetime

from cansnpscorer import vcf
from cansnpscorer.helpers.data_utils import (
    canSNP_extractor,
    nearest_distances,
    snp_extractor,
    variable_columns,
)
from cansnpscorer.helpers.file_utils import (
    read_alignment,
    read_genome_map,
    read_hierarchy,
    reference_parser,
)
from cansnpscorer.helpers.names import resolve_clusters
from cansnpscorer.hierarchy import HierarchyError
from cansnpscorer.model import AlignmentMap, ReferenceError
from cansnpscorer.validate import (
    check_alignment_against_reference,
    check_groups,
    report,
    unresolved_genomes,
)


def main(alignment=None, reference=None, hierarchy=None, genome_map=None, output=None,
         min_called=1.0) -> int:
    """Extract canSNPs from an alignment. Returns a process exit status."""
    logging.info(f"Started analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info(
        f"Parameters: alignment={alignment}, reference={reference}, hierarchy={hierarchy}, "
        f"genome_map={genome_map}, output={output}, min_called={min_called}")

    missing = [path for path in (alignment, reference, hierarchy, genome_map)
               if not os.path.isfile(path)]
    if missing:
        logging.error(f"These input files do not exist: {', '.join(missing)}")
        return 1

    try:
        tree = read_hierarchy(hierarchy)
        clusters = read_genome_map(genome_map)
        aligned = read_alignment(alignment)
        genome = reference_parser(reference)
    except (OSError, ValueError, HierarchyError, ReferenceError) as error:
        logging.error(error)
        return 1

    # Report everything inconsistent between the files before doing any work.
    problems, warnings = check_alignment_against_reference(aligned, genome)
    clusters, unmatched, ambiguous = resolve_clusters(clusters, list(aligned))
    found = unresolved_genomes(unmatched, ambiguous)
    problems, warnings = problems + found[0], warnings + found[1]
    found = check_groups(aligned, clusters, tree, next(iter(aligned)))
    problems, warnings = problems + found[0], warnings + found[1]
    if not report(problems, warnings):
        return 2

    try:
        isolates = tree.isolates(clusters)
        # The aligned's genome row carries the aligner's gaps, so the
        # column of a SNP is not its position in the genome genome.
        mapping = AlignmentMap(next(iter(aligned.values())), genome)
        snps = snp_extractor(aligned, min_called)
        # The flank of an assay is spoiled by any difference among the genomes
        # it will type, not only by canonical SNPs, so every column where the
        # scheme's genomes differ is found, and the outgroup's are left out.
        placed = sorted({name for names in clusters.values() for name in names})
        flank_columns = variable_columns(aligned, placed)
        cansnps = canSNP_extractor(aligned, isolates, snps, min_called)
    except (ValueError, HierarchyError, ReferenceError) as error:
        logging.error(error)
        return 1

    records: list[vcf.Record] = []
    gapped = 0
    for cansnp in cansnps:
        placed = mapping.locate(cansnp.column)
        if placed is None:
            # The reference is gapped in this column, so the SNP has no
            # position in the reference genome and no assay can target it.
            gapped += 1
            continue
        contig, position = placed
        if cansnp.is_canonical:
            records.append(vcf.snp_record(genome, contig.name, position,
                                          (cansnp.derived, cansnp.ancestral),
                                          cansnp.node, cansnp.ancestral))
        else:
            records.append(vcf.snp_record(genome, contig.name, position, cansnp.bases))
        records[-1].info["MISSING"] = str(cansnp.missing)

    if gapped:
        logging.warning(
            f"{gapped} SNPs sit in columns where the reference is gapped and were left out")

    variable: dict[str, list[int]] = {}
    for column in flank_columns:
        placed_at = mapping.locate(column)
        if placed_at is not None:
            variable.setdefault(placed_at[0].name, []).append(placed_at[1] - 1)
    for found in variable.values():
        found.sort()
    distances = nearest_distances([(record.chrom, record.pos - 1) for record in records],
                                  variable)
    for record, distance in zip(records, distances, strict=True):
        record.info["NEAREST"] = str(min(distance, vcf.NEAREST_CAP)) if distance is not None \
            else str(vcf.NEAREST_CAP)
    logging.info(f"{len(flank_columns)} columns differ among the {len(placed)} genomes of the "
                 f"scheme, which gives each SNP the distance to its nearest variation")

    header = vcf.header_lines(genome, reference, "extract",
                              [vcf.NODE, vcf.ANCESTRAL, vcf.NEAREST, vcf.MISSING])
    written = vcf.write(output, header, records)

    found = sum(1 for cansnp in cansnps if cansnp.is_canonical)
    without = sorted(set(tree.nodes) - {cansnp.node for cansnp in cansnps if cansnp.node})
    if without:
        logging.warning(
            f"{len(without)} of {len(tree)} nodes have no canSNP: {', '.join(without)}")
    logging.info(f"Wrote {written} SNPs to {output}, of which {found} are canonical")
    if min_called < 1.0:
        partial = sum(1 for cansnp in cansnps if cansnp.is_canonical and cansnp.missing)
        logging.info(f"{partial} of the canonical SNPs lie in columns where some genomes have "
                     f"no base (at least {min_called:g} of the genomes called); INFO/MISSING "
                     f"gives the number")

    logging.info(f"Finished analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    return 0
