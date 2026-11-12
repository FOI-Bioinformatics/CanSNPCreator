"""Extraction of canSNPs from a whole genome alignment."""

import logging
import os
from datetime import datetime
from typing import List

from cansnpscorer import vcf
from cansnpscorer.helpers.data_utils import canSNP_extractor, snp_extractor
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


def main(alignment=None, reference=None, hierarchy=None, genome_map=None, output=None) -> int:
    """Extract canSNPs from an alignment. Returns a process exit status."""
    logging.info(f"Started analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info(
        f"Parameters: alignment={alignment}, reference={reference}, hierarchy={hierarchy}, "
        f"genome_map={genome_map}, output={output}")

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
    problems += unresolved_genomes(unmatched, ambiguous)
    found = check_groups(aligned, clusters, tree, next(iter(aligned)))
    problems, warnings = problems + found[0], warnings + found[1]
    if not report(problems, warnings):
        return 2

    try:
        isolates = tree.isolates(clusters)
        # The aligned's genome row carries the aligner's gaps, so the
        # column of a SNP is not its position in the genome genome.
        mapping = AlignmentMap(next(iter(aligned.values())), genome)
        snps = snp_extractor(aligned)
        cansnps = canSNP_extractor(aligned, isolates, snps)
    except (ValueError, HierarchyError, ReferenceError) as error:
        logging.error(error)
        return 1

    records: List[vcf.Record] = []
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

    if gapped:
        logging.warning(
            f"{gapped} SNPs sit in columns where the reference is gapped and were left out")

    header = vcf.header_lines(genome, reference, "extract", [vcf.NODE, vcf.ANCESTRAL])
    written = vcf.write(output, header, records)

    found = sum(1 for cansnp in cansnps if cansnp.is_canonical)
    without = sorted(set(tree.nodes) - {cansnp.node for cansnp in cansnps if cansnp.node})
    if without:
        logging.warning(
            f"{len(without)} of {len(tree)} nodes have no canSNP: {', '.join(without)}")
    logging.info(f"Wrote {written} SNPs to {output}, of which {found} are canonical")

    logging.info(f"Finished analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    return 0
