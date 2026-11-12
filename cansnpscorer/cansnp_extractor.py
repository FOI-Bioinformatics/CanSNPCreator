"""Extraction of canSNPs from a whole genome alignment."""

import logging
import os
from datetime import datetime
from typing import List

from cansnpscorer.helpers.data_utils import canSNP_extractor, snp_extractor
from cansnpscorer.helpers.file_utils import (
    HCParser,
    cansnp_file_writer,
    read_alignment,
    read_hc_file,
    recognize_group_file_format,
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


def main(groupfile=None, hcfile=None, fastafile=None, referencefile=None, output=None) -> int:
    """Extract canSNPs from an alignment. Returns a process exit status."""
    logging.info(f"Started analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info(
        f"Parameters: groupfile={groupfile}, hcfile={hcfile}, fastafile={fastafile}, "
        f"referencefile={referencefile}, output={output}")

    missing = [path for path in (groupfile, hcfile, fastafile, referencefile)
               if not os.path.isfile(path)]
    if missing:
        logging.error(f"These input files do not exist: {', '.join(missing)}")
        return 1

    try:
        hierarchy = read_hc_file(hcfile)
        clusters = HCParser(groupfile, data_format=recognize_group_file_format(groupfile)).isolates
        alignment = read_alignment(fastafile)
        reference = reference_parser(referencefile)
    except (OSError, ValueError, HierarchyError, ReferenceError) as error:
        logging.error(error)
        return 1

    # Report everything inconsistent between the files before doing any work.
    problems, warnings = check_alignment_against_reference(alignment, reference)
    clusters, unmatched, ambiguous = resolve_clusters(clusters, list(alignment))
    problems += unresolved_genomes(unmatched, ambiguous)
    found = check_groups(alignment, clusters, hierarchy, next(iter(alignment)))
    problems, warnings = problems + found[0], warnings + found[1]
    if not report(problems, warnings):
        return 2

    try:
        isolates = hierarchy.isolates(clusters)
        # The alignment's reference row carries the aligner's gaps, so the
        # column of a SNP is not its position in the reference genome.
        mapping = AlignmentMap(next(iter(alignment.values())), reference)
        snps = snp_extractor(alignment)
        cansnps = canSNP_extractor(alignment, isolates, snps)
    except (ValueError, HierarchyError, ReferenceError) as error:
        logging.error(error)
        return 1

    rows: List[List] = []
    gapped = 0
    for cansnp in cansnps:
        placed = mapping.locate(cansnp.column)
        if placed is None:
            # The reference is gapped in this column, so the SNP has no
            # position in the reference genome and no assay can target it.
            gapped += 1
            continue
        contig, position = placed
        bases = (cansnp.derived, cansnp.ancestral) if cansnp.is_canonical else cansnp.bases
        rows.append([contig.name, position, cansnp.node or "None", *bases])

    if gapped:
        logging.warning(
            f"{gapped} SNPs sit in columns where the reference is gapped and were left out")

    cansnp_file_writer(rows, output)

    found = sum(1 for cansnp in cansnps if cansnp.is_canonical)
    without = sorted(set(hierarchy.nodes) - {cansnp.node for cansnp in cansnps if cansnp.node})
    if without:
        logging.warning(
            f"{len(without)} of {len(hierarchy)} nodes have no canSNP: {', '.join(without)}")
    logging.info(f"Wrote {len(rows)} SNPs, of which {found} are canonical")

    logging.info(f"Finished analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    return 0
