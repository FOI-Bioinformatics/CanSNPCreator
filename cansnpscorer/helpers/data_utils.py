import heapq
import logging
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple, Union

from cansnpscorer.helpers.genetic_code import codon_table

#: CanSNPer2 calls a node only when at least this many SNPs on its path are
#: seen, which is its --min_required_hits default.
MINIMUM_SNPS_PER_NODE = 3


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
    if not 0 < goal < 1:
        raise ValueError(f"the GC goal must lie between 0 and 1, got {goal}")

    # Extract the sequence around the SNP, clamped to the ends of the sequence.
    # Slicing with a negative start would otherwise wrap and give a window that
    # is short, in the wrong place, or empty.
    sequence = reference[max(snp - distance, 0):snp + distance]
    if not sequence:
        return 0.0

    gc_count = sequence.count('G') + sequence.count('C')
    observed = gc_count / float(len(sequence))

    # Normalized deviation from the target, clamped so that a sequence far from
    # a low goal cannot contribute a negative score.
    normalized_deviation = 1 - (abs(observed - goal) / goal)

    return max(normalized_deviation, 0.0)

#Binary search algorithm
def search(snp: int, intervals: Sequence[Sequence[int]], limit: int) -> Tuple[float, Optional[Sequence[int]]]:
    """
    Find the interval that contains the given position, or lies within limit of it.

    Args:
        snp (int): The position to search for.
        intervals (Sequence[Sequence[int]]): Intervals sorted by their start, each
            beginning with a start and an end. Further elements, such as a strand,
            are carried through untouched.
        limit (int): How far outside an interval still counts as close to it.

    Returns:
        Tuple[float, Optional[Sequence[int]]]: The normalized distance to the
        nearest interval and that interval, or (1.0, None) if none is close.
        The distance is 0.0 inside an interval and 1.0 at or beyond the limit.
        The second element is the interval only when the position is inside it,
        so that the caller knows the reading frame applies.
    """
    lowerbound = 0
    upperbound = len(intervals) - 1

    while lowerbound <= upperbound:
        current = (lowerbound + upperbound) // 2
        temp = intervals[current]

        left = max(temp[0] - limit, 0)
        right = temp[1] + limit

        if left <= snp <= right:
            if snp < temp[0]:
                return (temp[0] - snp) / float(limit), None
            if snp > temp[1]:
                return (snp - temp[1]) / float(limit), None
            return 0.0, temp
        if snp < left:
            upperbound = current - 1
        else:
            lowerbound = current + 1

    return 1.0, None

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

@dataclass
class CanSnp:
    """One variable column of the alignment, and the node it marks, if any.

    ``derived`` is the base carried by every isolate below the node and by no
    other genome, which is what makes the SNP canonical; ``ancestral`` is the
    other base. They are only known for a canonical SNP, so both are None when
    ``node`` is.
    """

    column: int
    bases: Tuple[str, ...]
    node: Optional[str] = None
    derived: Optional[str] = None
    ancestral: Optional[str] = None
    also_matched: Tuple[str, ...] = ()
    """Further nodes separated by the same column, which makes it ambiguous."""

    @property
    def is_canonical(self) -> bool:
        return self.node is not None


def canSNP_extractor(alignment: Dict[str, str], groups: Dict[str, List[str]],
                     snps: List[Dict[str, Union[int, Dict[str, int]]]]) -> List[CanSnp]:
    """
    Find which of the SNPs separate exactly the isolates below one node.

    A SNP is canonical for a node when every isolate below that node carries
    the same base and no genome outside it does. The base they carry is the
    derived one, so the polarity of the SNP follows from the node it marks.

    Args:
        alignment: The aligned sequences, keyed by name.
        groups: The isolates below each node, which for an internal node means
            the isolates of every leaf beneath it.
        snps: The variable columns, as snp_extractor reports them.

    Returns:
        One record per SNP, in the order given, canonical or not.

    Raises:
        ValueError: If a group names an isolate that the alignment does not hold,
            which would otherwise be a KeyError part way through the run.
    """
    logging.info(f"Looking for canSNPs among {len(snps)} SNPs over {len(groups)} nodes")

    wanted = {isolate for isolates in groups.values() for isolate in isolates}
    missing = sorted(wanted - set(alignment))
    if missing:
        raise ValueError(
            f"the group file names {len(missing)} isolates that the alignment does not hold: "
            f"{', '.join(missing[:5])}")

    # Nodes with no isolates cannot be marked by any SNP. They are reported
    # once here rather than silently never matching.
    empty = sorted(node for node, isolates in groups.items() if not isolates)
    if empty:
        logging.warning(
            f"{len(empty)} nodes have no isolates, so no canSNP can mark them: "
            f"{', '.join(empty[:5])}")

    members = {
        node: list(isolates)
        for node, isolates in groups.items() if isolates
    }

    cansnps: List[CanSnp] = []
    canonical = 0
    ambiguous = 0
    complementary = 0

    for snp in snps:
        column = snp['index']
        counts = snp['variants']
        bases = tuple(counts)

        matches: List[Tuple[str, str]] = []
        for node, isolates in members.items():
            shared = alignment[isolates[0]][column]
            if any(alignment[isolate][column] != shared for isolate in isolates[1:]):
                continue
            if counts.get(shared, 0) == len(isolates):
                matches.append((node, shared))

        if not matches:
            cansnps.append(CanSnp(column=column, bases=bases))
            continue

        node, derived = matches[0]
        other = [base for base in bases if base != derived]
        cansnps.append(
            CanSnp(
                column=column,
                bases=bases,
                node=node,
                derived=derived,
                ancestral=other[0] if other else None,
                also_matched=tuple(name for name, _base in matches[1:]),
            )
        )
        canonical += 1
        if len(matches) > 1:
            first = set(members[matches[0][0]])
            if any(set(members[other]) == first for other, _base in matches[1:]):
                ambiguous += 1
            else:
                complementary += 1

    if ambiguous:
        logging.warning(
            f"{ambiguous} canSNPs separate more than one node, because those nodes hold the "
            "same isolates; each is reported against the first of them")
    if complementary:
        logging.warning(
            f"{complementary} canSNPs separate two nodes whose isolates together are every "
            f"genome of the alignment, so which side carries the derived base cannot be told "
            f"without a genome outside both; each is reported against the first of them. An "
            f"outgroup genome in the alignment, placed at no node, settles it")
    logging.info(f"Found {canonical} canSNPs")

    return cansnps


def apply_filters(scores, cutoff, top=None):
    """
    Mark the scores that fail the cutoff, or are not among the best of their node.

    Nothing is dropped: a score at or below the cutoff gains the FILTER code
    LowScore, and one that passes it but is not among the ``top`` best of its
    node gains NotTop, so that the output still shows what was considered.

    Args:
        scores: Scored SNPs, each with a total, a node and a list of filters.
        cutoff: A score must be above this to pass.
        top: How many of the best scores of each node pass. All of them if None.

    Returns:
        The same scores, in the order given.
    """
    for score in scores:
        if score.total <= cutoff and "LowScore" not in score.filters:
            score.filters.append("LowScore")

    if top is not None:
        by_node = defaultdict(list)
        for position, score in enumerate(scores):
            if not score.filters:
                # A min-heap, so the totals are negated to pop the best first; the
                # position breaks ties without comparing the scores themselves.
                heapq.heappush(by_node[score.node], (-score.total, position, score))
        for node_scores in by_node.values():
            for _ in range(min(len(node_scores), top)):
                heapq.heappop(node_scores)
            for _total, _position, score in node_scores:
                score.filters.append("NotTop")

    return scores


def summarize_scores(nodes_data, scores):
    """
    Count the SNPs kept for each node, and report the nodes left without one.

    Args:
        nodes_data: The node identifiers of the scheme.
        scores: The scored SNPs that passed filtering.

    Returns:
        dict: The nodes that have at least one SNP, and how many.
    """
    node_counts = {node: 0 for node in nodes_data}

    for score in scores:
        if score.node in node_counts:
            node_counts[score.node] += 1

    represented = {node: count for node, count in node_counts.items() if count}
    for node, count in represented.items():
        logging.info(f"{node}: {count} canSNPs")

    without = [node for node, count in node_counts.items() if not count]
    if without:
        logging.warning(
            f"{len(without)} of {len(node_counts)} nodes have no canSNP above the cutoff, so "
            f"they cannot be typed: {', '.join(without)}")

    # CanSNPer2 calls a node only when enough of the SNPs on its path are
    # observed, three by default, so a node with fewer is worth naming.
    shallow = [node for node, count in represented.items() if count < MINIMUM_SNPS_PER_NODE]
    if shallow:
        logging.warning(
            f"{len(shallow)} nodes have fewer than {MINIMUM_SNPS_PER_NODE} canSNPs, which is "
            f"what CanSNPer2 requires on a path by default: {', '.join(shallow)}")

    return represented
