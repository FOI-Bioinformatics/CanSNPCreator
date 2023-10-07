from typing import List, Dict, Union, NoReturn, Tuple
import re
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
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


def snp_extractor_worker(worker_id, start, end, alignment, reference, nucleotides):
    logging.info(f"Worker {worker_id} started, processing indices {start} to {end}")
    snps = []
    for index in range(start, end):
        nucleotide = reference[index]
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

    logging.info(f"Worker {worker_id} completed")
    return snps

def snp_extractor(fasta: Dict[str, str], num_threads=4) -> List[Dict[str, Union[int, Dict[str, int]]]]:
    alignment = list(fasta.values())
    reference = alignment[0]
    alignment = alignment[1:]

    nucleotides = set(['A', 'T', 'C', 'G'])
    total_length = len(reference)

    # Split the task into chunks for each thread
    chunk_size = total_length // num_threads
    futures = []

    logging.info(f"Starting SNP extraction with {num_threads} workers...")

    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        for i in range(num_threads):
            start_index = i * chunk_size
            end_index = (i + 1) * chunk_size if i != num_threads - 1 else total_length
            futures.append(executor.submit(snp_extractor_worker, i, start_index, end_index, alignment, reference, nucleotides))

    # Collect the results
    snps = []
    for future in futures:
        snps.extend(future.result())

    # Sort the SNPs by their index
    snps = sorted(snps, key=lambda x: x['index'])

    logging.info(f"Extraction of {len(snps)} SNPs complete")
    return snps


def canSNP_extractor_worker(worker_id, start, end, snps, alignment, groups):
    logging.info(f"Worker {worker_id} started, processing indices {start} to {end}")
    cansnps = []
    canSNP_count = 0  # Counter for the number of canSNPs

    for i in range(start, end):
        snp = snps[i]
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

    logging.info(f"Worker {worker_id} completed. Total canSNPs extracted in this worker: {canSNP_count}")

    return cansnps


def canSNP_extractor(alignment: Dict[str, str], groups: Dict[str, List[str]],
                     snps: List[Dict[str, Union[int, Dict[str, int]]]], num_threads=4) -> List[List[Union[str, int]]]:
    logging.info(f"Starting canSNP extraction with {num_threads} workers...")

    total_snps = len(snps)

    # Split the task into chunks for each thread
    chunk_size = total_snps // num_threads
    futures = []

    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        for i in range(num_threads):
            start_index = i * chunk_size
            end_index = (i + 1) * chunk_size if i != num_threads - 1 else total_snps
            futures.append(executor.submit(canSNP_extractor_worker, i, start_index, end_index, snps, alignment, groups))

    # Collect the results
    cansnps = []
    for future in futures:
        cansnps.extend(future.result())

    # Sort based on the SNP index, which is the first element in each sublist
    cansnps.sort(key=lambda x: int(x[0]))

    logging.info(f"CanSNP extraction complete. Total canSNPs extracted: {len(cansnps)}")

    return cansnps