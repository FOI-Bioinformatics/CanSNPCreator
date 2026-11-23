"""A sample of SNP sites from several alignments, as one matrix for a tree.

The genomes a scheme was made from and the genomes placed on it later are
often aligned separately against the same reference: RepGenR2's alignment of
the representatives and an alignment of further genomes. A tree over all of
them needs one matrix whose columns are the same reference positions in every
row. ``sample_sites`` takes a random sample of the SNP positions extract
found, maps each position to a column of every alignment with that
alignment's own ``AlignmentMap``, so that alignments with different gap
columns still agree, and keeps only those columns of each record.

The sample is that of Python's ``random.Random(seed).sample`` over the SNP
positions in the order of the VCF file, so the same file, number of sites and
seed give the same columns.
"""

import logging
import random
from collections.abc import Sequence
from itertools import chain

from cansnpscorer.helpers.file_utils import iter_fasta, open_text, reference_parser
from cansnpscorer.model import AlignmentMap, Reference, ReferenceError

#: The character written where an alignment has no column for a position.
NO_COLUMN = "-"


def snp_positions(snps_vcf: str, reference: Reference) -> list[int]:
    """The 0-based positions, over the concatenated reference, of the SNPs of
    a VCF file, in file order and each once.

    Only the CHROM and POS columns are read, so any VCF naming the reference's
    contigs will do, such as the one extract writes.

    Raises:
        ReferenceError: If a SNP names a contig the reference lacks, or lies
            outside its contig.
        ValueError: If a record line cannot be read.
    """
    positions: dict[int, None] = {}
    with open_text(snps_vcf) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t", 2)
            if len(fields) < 2:
                raise ValueError(f"{snps_vcf}: line {number} is not a VCF record")
            try:
                position = int(fields[1])
            except ValueError:
                raise ValueError(f"{snps_vcf}: line {number} gives the position '{fields[1]}'") from None
            contig = reference.contig(fields[0])
            if not 1 <= position <= contig.length:
                raise ReferenceError(
                    f"{snps_vcf}: line {number} gives position {position}, outside contig "
                    f"'{contig.name}' of {contig.length} bases"
                )
            positions.setdefault(contig.start + position - 1, None)
    return list(positions)


def sample_sites(
    snps_vcf: str, reference: Reference, alignments: Sequence[str], sites: int = 20000, seed: int = 0
) -> dict[str, str]:
    """Sample SNP sites and take them from every record of the alignments.

    Args:
        snps_vcf: The SNPs to sample from, as extract writes them.
        reference: The reference genome every alignment is against.
        alignments: Whole genome alignments, each with the reference as its
            first record, gapped or projected as ``AlignmentMap`` accepts.
        sites: How many SNP positions to sample; all of them when there are
            fewer.
        seed: The seed of the random sample.

    Returns:
        One row per record, its bases at the sampled positions in reference
        order. A position an alignment has no column for is a gap in that
        alignment's rows. A record in two alignments is taken from the first.

    Raises:
        ValueError: If an alignment's records differ in length, or it holds
            no record.
        ReferenceError: If an alignment's first record is not the reference.
    """
    positions = snp_positions(snps_vcf, reference)
    taken = sorted(random.Random(seed).sample(positions, min(sites, len(positions))))
    logging.info(f"Sampled {len(taken)} of {len(positions)} SNP positions with seed {seed}")

    rows: dict[str, str] = {}
    for path in alignments:
        records = iter_fasta(path)
        first = next(records, None)
        if first is None:
            raise ValueError(f"{path}: holds no FASTA records")
        reference_name, aligned_reference = first
        mapping = AlignmentMap(aligned_reference, reference)
        columns = [mapping.to_column(position) for position in taken]
        missing = sum(column is None for column in columns)
        if missing:
            logging.warning(
                f"{path}: {missing} of the sampled positions have no column, and are written as gaps"
            )
        width = len(aligned_reference)
        added = repeated = 0
        names: set[str] = set()
        # Chained rather than unpacked, so that one record is held at a time.
        for name, sequence in chain([(reference_name, aligned_reference)], records):
            if len(sequence) != width:
                raise ValueError(
                    f"{path}: the records are not aligned; '{name}' holds {len(sequence)} "
                    f"columns and the reference row {width}"
                )
            if name in names:
                raise ValueError(
                    f"{path}: the name '{name}' is used more than once; sequence IDs, the "
                    "first word of each header, must be unique"
                )
            names.add(name)
            if name in rows:
                repeated += 1
                continue
            rows[name] = "".join(NO_COLUMN if column is None else sequence[column] for column in columns)
            added += 1
        logging.info(
            f"Took {added} records from {path}"
            + (f"; {repeated} were taken from an earlier alignment" if repeated else "")
        )
    return rows


def main(snps=None, reference=None, alignments=(), sites=20000, seed=0, output=None) -> int:
    """Write sampled SNP sites of the alignments as FASTA. Returns an exit status."""
    try:
        rows = sample_sites(snps, reference_parser(reference), alignments, sites, seed)
        with open(output, "w", newline="\n") as handle:
            handle.writelines(f">{name}\n{row}\n" for name, row in rows.items())
    except (OSError, ValueError) as error:
        logging.error(error)
        return 1
    logging.info(f"Wrote {len(rows)} records to {output}")
    return 0
