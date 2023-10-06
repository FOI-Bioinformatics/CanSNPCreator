from typing import List, Dict, Union, NoReturn, Tuple
import re
from collections import Counter
import logging
from cansnpscorer.helpers.genetic_code import codon_table

def neighbour_distance(snp: int, left: int, right: int, min_dist: int = 25) -> float:
    """
    Calculate the normalized distance of a SNP to its neighboring SNPs.

    Args:
        snp (int): The position of the SNP of interest.
        left (int): The position of the left neighboring SNP.
        right (int): The position of the right neighboring SNP.
        min_dist (int, optional): The minimum distance to be considered. Defaults to 25.

    Returns:
        float: The normalized distance to the closest neighbor.
    """
    left_dist = min(abs(left - snp), min_dist)
    right_dist = min(abs(right - snp), min_dist)

    max_dist = min(left_dist, right_dist)

    return max_dist / float(min_dist)

#Calculate gc content
def gc_content(snp: int, reference: str, goal: float, distance: int = 50) -> float:
    """
    Calculate the GC content around a given SNP and normalize it relative to a goal GC content.

    Args:
        snp (int): The position of the SNP in the reference.
        reference (str): The reference genome sequence.
        goal (float): The target GC content.
        distance (int, optional): The distance to look on either side of the SNP. Defaults to 50.

    Returns:
        float: Normalized deviation from the target GC content.
    """
    # Extract the sequence around the SNP
    sequence = reference[snp - distance:snp + distance]

    # Count the number of G's and C's in the sequence
    gc_count = sequence.count('G') + sequence.count('C')

    # Calculate the GC content
    gc_content = gc_count / float(len(sequence))

    # Calculate the normalized deviation from the target GC content
    normalized_deviation = 1 - (abs(gc_content - goal) * (1 / goal))

    return normalized_deviation

#Binary search algorithm
def search(snp: int, intervals: list[Tuple[int, int]], limit: int) -> Union[Tuple[int], Tuple[float, Tuple[int, int]]]:
    """
    Perform a binary search to find the interval containing the given SNP.

    Args:
        snp (int): The position of the SNP to search for.
        intervals (list[Tuple[int, int]]): A list of sorted intervals to search within.
        limit (int): A limit used for extending the interval boundaries.

    Returns:
        Union[Tuple[int], Tuple[float, Tuple[int, int]]]: Either a tuple containing just 1 if not found,
        or a tuple containing the normalized distance and the interval if found.
    """
    lowerbound = 0
    upperbound = len(intervals) - 1

    while True:
        if upperbound < lowerbound:
            return 1,

        current = (lowerbound + upperbound) // 2  # Use // for integer division in Python 3
        temp = intervals[current]

        left = max(temp[0] - limit, 0)
        right = temp[1] + limit

        if left <= snp <= right:
            if snp < temp[0]:
                return (temp[0] - snp) / float(limit),
            elif snp > temp[1]:
                return (snp - temp[1]) / float(limit),
            else:
                return 0, temp
        elif snp < left:
            upperbound = current - 1
        else:  # snp > right
            lowerbound = current + 1

#Method for checking type of substitution mutation
def gene_checker(snp: List[Union[int, str]], reference: List[str], gene_interval: Tuple[int, int, str], correction: int) -> float:
    """
    Checks the gene region for SNPs and compares the codons.

    Args:
        snp (List[Union[int, str]]): List representing the SNP, where the first element is the index and the rest are variants.
        reference (List[str]): List representing the reference genome sequence.
        gene_interval (Tuple[int, int, str]): Tuple representing the start, end, and strand of the gene.
        correction (int): Position correction factor.

    Returns:
        float: Similarity score between the codons.
    """
    variants = []
    index = snp[0] - gene_interval[0] - correction
    gene = reference[0][gene_interval[0] + correction: gene_interval[1] + correction + 1]

    if gene_interval[2] == "-":
        index = (gene_interval[1] - gene_interval[0]) - index
        gene = gene[::-1]

    j = index % 3
    for i in snp[2:]:
        left = gene[index - j: index]
        mid = i
        right = gene[index + 1: index + 1 + 2 - j]

        variants.append(left + mid + right)

    sim_score = codon_comparer(variants, gene_interval[2])
    return sim_score


#Compare two codons
def codon_comparer(variants: List[str], strand: str) -> float:
    """
    Compare two codons to determine their similarity.

    Args:
        variants (List[str]): List containing the two codons to compare.
        strand (str): DNA strand orientation ('+' or '-').

    Returns:
        float: Similarity score between the codons.
    """
    codon1, codon2 = variants[0], variants[1]

    if strand == "-":
        codon1 = codon_complement(codon1)
        codon2 = codon_complement(codon2)

    amino_acid1 = codon_table.get(codon1, '')
    amino_acid2 = codon_table.get(codon2, '')

    if amino_acid1 == amino_acid2:
        return 1.0
    elif amino_acid1 == "*" or amino_acid2 == "*":
        return 0.0
    else:
        return 0.5


#Translates codon in to their complement
def codon_complement(codon: str) -> str:
    """
    Translates a codon into its complement.

    Args:
        codon (str): The original codon.

    Returns:
        str: The complement codon.
    """
    translation = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join([translation[nucleotide] for nucleotide in codon])


def clean_name(name: str) -> str:
    """
    Cleans the input name by removing all non-alphanumeric characters.

    Args:
        name (str): The original name string that may contain non-alphanumeric characters.

    Returns:
        str: The cleaned name containing only alphanumeric characters.
    """
    return ''.join(filter(str.isalnum, name))


def snp_extractor(fasta: Dict[str, str]) -> List[Dict[str, Union[int, Dict[str, int]]]]:
    """
    Extract Single Nucleotide Polymorphisms (SNPs) from genomic alignments.

    Args:
        fasta (Dict[str, str]): A dictionary where keys are sequence names and values are genomic sequences.

    Returns:
        List[Dict[str, Union[int, Dict[str, int]]]]: A list of dictionaries. Each dictionary contains:
            - 'index': The position of the SNP in the reference sequence.
            - 'variants': A dictionary with the nucleotide variants at that position and their counts.

    Example:
        Return value might look like:
        [
            {'index': 5, 'variants': {'A': 3, 'T': 2}},
            {'index': 9, 'variants': {'G': 2, 'C': 3}},
            ...
        ]

    """
    alignment = list(fasta.values())
    reference = alignment[0]
    alignment = alignment[1:]
    snps = []

    nucleotides = set(['A', 'T', 'C', 'G'])
    total_length = len(reference)
    last_logged_percent = 0

    logging.info("Starting SNP extraction...")

    for index, nucleotide in enumerate(reference):
        if nucleotide not in nucleotides:
            continue

        variants_at_position = [genome[index] for genome in alignment]
        variants_count = Counter(variants_at_position)

        if nucleotide not in variants_count:
            variants_count[nucleotide] = 1
        else:
            variants_count[nucleotide] += 1

        if all(base in nucleotides for base in variants_count) and len(variants_count) == 2:
            snps.append({'index': index, 'variants': dict(variants_count)})

        # Log progress every 10%
        current_percent = int(((index + 1) / total_length) * 100)
        if current_percent >= last_logged_percent + 10:
            logging.info(f"Processed {index + 1}/{total_length} positions ({current_percent}%)")
            last_logged_percent = current_percent

    logging.info(f"Extraction of {len(snps)} SNPs complete")
    return snps

def canSNP_extractor(alignment: Dict[str, str], groups: Dict[str, List[str]], snps: List[Dict[str, Union[int, Dict[str, int]]]]) -> List[List[Union[str, int]]]:
    """
    Extract Canonical Single Nucleotide Polymorphisms (canSNPs) from genomic alignments.

    Args:
        alignment (Dict[str, str]): A dictionary where keys are sequence names and values are genomic sequences.
        groups (Dict[str, List[str]]): A dictionary defining the groups of sequences.
        snps (List[Dict[str, Union[int, Dict[str, int]]]]): A list of dictionaries, each containing 'index' and 'variants',
                                                             representing the SNPs extracted from the genomic alignments.

    Returns:
        List[List[Union[str, int]]]: A list containing the canonical SNPs. Each element is a list containing the SNP index,
                                     the group it belongs to, and the variants present at that position.

    Example:
        Return value might look like:
        [
            ['5', 'Group1', 'A', 'T'],
            ['9', 'Group2', 'G', 'C'],
            ...
        ]

    """
    logging.info("Starting canSNP extraction...")
    total_snps = len(snps)
    cansnps = []

    canSNP_count = 0  # Counter for the number of canSNPs
    last_logged_percent = -1  # Initialize to an impossible value

    for i, snp in enumerate(snps):
        variant = snp['index']
        frequencies = snp['variants']
        cansnp = "None"

        for group in groups:
            nucleotide = ""
            frequency = 1

            for index, isolate in enumerate(groups[group]):
                current = alignment[clean_name(isolate)][variant]

                if index == 0:
                    nucleotide = current
                else:
                    if nucleotide == current:
                        frequency += 1
                    else:
                        break
            else:
                if frequency == frequencies.get(nucleotide, 0):
                    cansnp = group

        cansnps.append([str(variant), cansnp] + list(frequencies.keys()))

        # Increment canSNP counter if cansnp is not "None"
        if cansnp != "None":
            canSNP_count += 1

        # Log progress every 10%
        current_percent = int(((i + 1) / total_snps) * 100)
        if current_percent >= last_logged_percent + 10:
            logging.info(f"Writing progress: {current_percent}%")
            last_logged_percent = current_percent  # Update the last logged percentage only when it changes

    logging.info(f"CanSNP extraction complete. Total canSNPs extracted: {canSNP_count}")

    return cansnps
