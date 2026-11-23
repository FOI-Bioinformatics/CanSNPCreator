"""How many times the sequence around a SNP occurs in the reference.

A SNP whose flank is repeated elsewhere in the genome, as in rRNA operons,
insertion sequences or paralogous genes, is harder to align to the right copy
in a query genome. On 1,229 Brucella genomes left out of a scheme's tree, such
canSNPs lacked a base three times as often as others, and disagreed with the
other canSNPs of their node three times as often.

For each SNP, the k-mers that hold it are looked up on both strands of every
contig, and the copy number of the SNP is the highest count among them, so a
repeat of at least k bases that spans the SNP is found. With a step above 1
only every step-th of those k-mers is looked up, which uses less memory and
finds every such repeat of at least k + step - 1 bases. K-mers holding a base
other than A, C, G or T, such as the N of a scaffold gap, are not counted.

The reference base is used at the SNP, so the count is that of the reference
allele's context.
"""

from collections.abc import Iterable

from cansnpscorer.model import Reference

BASES = frozenset("ACGT")
COMPLEMENT = str.maketrans("ACGT", "TGCA")

#: Default k-mer size: long enough to be unique by chance in a bacterial genome.
DEFAULT_K = 31
#: Default distance between the starts of the k-mers looked up around one SNP:
#: every k-mer that holds it.
DEFAULT_STEP = 1


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def copy_numbers(
    reference: Reference, sites: Iterable[tuple[str, int]], k: int = DEFAULT_K, step: int = DEFAULT_STEP
) -> dict[tuple[str, int], int]:
    """The copy number of the context of each site.

    Args:
        reference: The reference genome.
        sites: (contig name, 0-based position) of each SNP.
        k: The k-mer size.
        step: The distance between the starts of the k-mers looked up.

    Returns:
        For each site, the highest number of times one of its k-mers occurs on
        either strand of the reference; 1 when all are unique, and 0 when no
        k-mer of only A, C, G and T fits around the site.
    """
    if k < 1 or step < 1:
        raise ValueError(f"the k-mer size and step must be at least 1, got {k} and {step}")

    # Each k-mer looked up is keyed by the hash of its sequence; the hash of
    # its reverse complement points to the same entry, so that an occurrence
    # on either strand is counted while the forward strand alone is read. A
    # k-mer that is its own reverse complement has one key, and each of its
    # occurrences is counted once, as it should be.
    wanted: dict[int, int] = {}  # hash of a k-mer or its reverse complement -> entry
    entries_of: dict[tuple[str, int], tuple[int, ...]] = {}
    entries = 0
    for name, position in sites:
        sequence = reference.contig(name).sequence
        # Starts from k-1 bases before the SNP to the SNP itself, at most
        # step apart, the last one always at the SNP.
        starts = list(range(position - k + 1, position + 1, step))
        if starts[-1] != position:
            starts.append(position)
        found: list[int] = []
        for start in starts:
            if start < 0 or start + k > len(sequence):
                continue
            kmer = sequence[start : start + k].upper()
            if not set(kmer) <= BASES:
                continue
            key = hash(kmer)
            entry = wanted.get(key)
            if entry is None:
                entry = entries
                entries += 1
                wanted[key] = entry
                wanted[hash(_reverse_complement(kmer))] = entry
            found.append(entry)
        entries_of[(name, position)] = tuple(found)

    counts = [0] * entries
    if entries:
        for contig in reference.contigs:
            sequence = contig.sequence.upper()
            for start in range(len(sequence) - k + 1):
                entry = wanted.get(hash(sequence[start : start + k]))
                if entry is not None:
                    counts[entry] += 1

    return {site: max((counts[entry] for entry in found), default=0) for site, found in entries_of.items()}
