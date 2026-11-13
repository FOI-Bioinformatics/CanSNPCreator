import logging
import math
from bisect import bisect_right
from collections import Counter, defaultdict
from collections.abc import Sequence
from dataclasses import dataclass

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
def gc_content(snp: int, reference: str, goal: float, distance: int = 50,
               tolerance: float = 0.0) -> float:
    """
    Calculate the GC content around a given SNP and normalize it relative to a goal GC content.

    Args:
        snp (int): The position of the SNP in the reference.
        reference (str): The reference genome sequence.
        goal (float): The target GC content.
        distance (int, optional): The distance to look on either side of the SNP. Defaults to 50.
        tolerance (float, optional): A deviation from the goal up to this counts as none,
            as primers are designed within a range of GC content. Defaults to 0.

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
    deviation = max(abs(observed - goal) - tolerance, 0.0)
    normalized_deviation = 1 - (deviation / goal)

    return max(normalized_deviation, 0.0)

#Binary search algorithm
def search(snp: int, intervals: Sequence[Sequence[int]], limit: int) -> tuple[float, Sequence[int] | None]:
    """
    Find the interval that holds a position, or the nearest one within limit of it.

    Args:
        snp: The position.
        intervals: Intervals sorted by their start that do not overlap, each
            beginning with a start and an end, inclusive. Further elements, such
            as a strand, are carried through.
        limit: How far outside an interval still counts as near it.

    Returns:
        The distance to the nearest interval as a fraction of limit, 0.0 inside
        one and 1.0 at or beyond the limit, and that interval when the position
        is inside it, or None.
    """
    index = bisect_right([interval[0] for interval in intervals], snp) - 1
    nearest = None
    # Only the interval starting at or before the position, and the one after
    # it, can be the nearest when the intervals do not overlap.
    for candidate in (index, index + 1):
        if not 0 <= candidate < len(intervals):
            continue
        start, end = intervals[candidate][0], intervals[candidate][1]
        if start <= snp <= end:
            return 0.0, intervals[candidate]
        distance = start - snp if snp < start else snp - end
        if nearest is None or distance < nearest:
            nearest = distance
    if nearest is None or limit <= 0 or nearest >= limit:
        return 1.0, None
    return nearest / float(limit), None


def codon_comparer(variants: list[str], strand: str) -> float:
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
    translation = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join([translation[nucleotide] for nucleotide in codon])


def called_needed(total: int, min_called: float) -> int:
    """The number of genomes out of ``total`` that ``min_called`` asks for."""
    return max(1, math.ceil(min_called * total - 1e-9))


def snp_extractor(fasta: dict[str, str],
                  min_called: float = 1.0) -> list[dict[str, int | dict[str, int]]]:
    """
    Extract Single Nucleotide Polymorphisms (SNPs) from genomic alignments.

    Args:
        fasta (Dict[str, str]): A dictionary where keys are sequence names and values are genomic sequences.
        min_called: The fraction of genomes that must carry a base (A, C, G or
            T) in a column for it to count. With 1.0 a single gap or N in any
            genome leaves the column out; below it the column is a SNP when
            the genomes that are called carry exactly two bases.

    Returns:
        List[Dict[str, Union[int, Dict[str, int]]]]: A list of dictionaries. Each dictionary contains:
            - 'index': The position of the SNP in the reference sequence.
            - 'variants': A dictionary with the nucleotide variants at that position and their counts.
            - 'missing': The number of genomes without a base in the column.

    Example:
        Return value might look like:
        [
            {'index': 5, 'variants': {'A': 3, 'T': 2}},
            {'index': 9, 'variants': {'G': 2, 'C': 3}},
            ...
        ]

    """
    nucleotides = frozenset("ACGT")
    rows = list(fasta.values())
    total_length = len(rows[0]) if rows else 0
    snps = []
    step = max(total_length // 10, 1)

    logging.info("Starting SNP extraction...")

    # zip walks the columns in C, one tuple per column, which is several times
    # faster than indexing every genome in Python, and the counts include the
    # reference row as the first element.
    needed = called_needed(len(rows), min_called)
    for index, column in enumerate(zip(*rows, strict=True)):
        if column[0] not in nucleotides:
            continue
        counts = Counter(column)
        if len(counts) == 2 and nucleotides.issuperset(counts):
            snps.append({'index': index, 'variants': dict(counts), 'missing': 0})
        elif needed < len(rows):
            called = {base: n for base, n in counts.items() if base in nucleotides}
            total = sum(called.values())
            if len(called) == 2 and total >= needed:
                snps.append({'index': index, 'variants': called,
                             'missing': len(rows) - total})
        if (index + 1) % step == 0:
            logging.info(f"Processed {index + 1}/{total_length} positions "
                         f"({round(100 * (index + 1) / total_length)}%)")

    logging.info(f"Extraction of {len(snps)} SNPs complete")
    return snps

def variable_columns(fasta: dict[str, str], genomes: Sequence[str]) -> list[int]:
    """
    The columns where at least two of these genomes carry different bases.

    Unlike snp_extractor this keeps columns where a third base occurs or where
    other genomes have no base, since any of them spoils the flank of an assay
    as much as a canonical SNP does. A gap is not counted as a difference,
    because in an alignment it cannot be told apart from missing sequence.

    Args:
        fasta: The aligned sequences, keyed by name.
        genomes: The genomes to compare, such as those placed in the scheme.
    """
    from operator import itemgetter

    names = list(fasta)
    indexes = [names.index(genome) for genome in genomes if genome in fasta]
    if len(indexes) < 2:
        return []
    pick = itemgetter(*indexes)
    nucleotides = frozenset("ACGT")
    found: list[int] = []
    for index, column in enumerate(zip(*fasta.values(), strict=True)):
        # Most columns are the same in every genome; set() on the whole column
        # rules them out before the chosen genomes are looked at.
        if len(set(column)) < 2:
            continue
        if len(set(pick(column)) & nucleotides) >= 2:
            found.append(index)
    return found


def nearest_distances(positions: Sequence[tuple[str, int]],
                      variable: dict[str, list[int]]) -> list[int | None]:
    """For each (contig, 0-based position), the distance to the nearest other
    variable position on the same contig, or None when there is none."""
    distances: list[int | None] = []
    for contig, position in positions:
        found = variable.get(contig, [])
        index = bisect_right(found, position)
        candidates = [abs(found[i] - position) for i in (index - 2, index - 1, index)
                      if 0 <= i < len(found) and found[i] != position]
        distances.append(min(candidates) if candidates else None)
    return distances


@dataclass
class CanSnp:
    """One variable column of the alignment, and the node it marks, if any.

    ``derived`` is the base carried by every isolate below the node and by no
    other genome, which is what makes the SNP canonical; ``ancestral`` is the
    other base. They are only known for a canonical SNP, so both are None when
    ``node`` is.
    """

    column: int
    bases: tuple[str, ...]
    node: str | None = None
    derived: str | None = None
    ancestral: str | None = None
    also_matched: tuple[str, ...] = ()
    """Further nodes separated by the same column, which makes it ambiguous."""
    missing: int = 0
    """Genomes without a base (gap or N) in the column."""

    @property
    def is_canonical(self) -> bool:
        return self.node is not None


def canSNP_extractor(alignment: dict[str, str], groups: dict[str, list[str]],
                     snps: list[dict[str, int | dict[str, int]]],
                     min_called: float = 1.0) -> list[CanSnp]:
    """
    Find which of the SNPs separate exactly the isolates below one node.

    A SNP is canonical for a node when every isolate below that node carries
    the same base and no genome outside it does. The base they carry is the
    derived one, so the polarity of the SNP follows from the node it marks.

    In a column with missing bases (``min_called`` below 1.0) only the genomes
    with a base are compared: at least ``min_called`` of the node's isolates
    must be called, all of them with the derived base, and no called genome
    outside the node may carry it.

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
            f"the genome map names {len(missing)} isolates that the alignment does not hold: "
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

    cansnps: list[CanSnp] = []
    canonical = 0
    ambiguous = 0
    complementary = 0
    undecided = 0

    for snp in snps:
        column = snp['index']
        counts = snp['variants']
        bases = tuple(counts)

        missing = snp.get('missing', 0)
        matches: list[tuple[str, str]] = []
        for node, isolates in members.items():
            if not missing:
                shared = alignment[isolates[0]][column]
                if any(alignment[isolate][column] != shared for isolate in isolates[1:]):
                    continue
                if counts.get(shared, 0) == len(isolates):
                    matches.append((node, shared))
                continue
            called = [alignment[isolate][column] for isolate in isolates]
            called = [base for base in called if base in counts]
            if len(called) < called_needed(len(isolates), min_called):
                continue
            shared = called[0]
            if any(base != shared for base in called[1:]):
                continue
            if counts.get(shared, 0) == len(called):
                matches.append((node, shared))

        if len(matches) > 1 and missing:
            # A node and one below it can both match only when the genomes
            # that tell them apart have no base here, so the column does not
            # say which of the two it marks.
            sets = [set(members[name]) for name, _base in matches]
            if any(a < b for a in sets for b in sets):
                undecided += 1
                matches = []

        if not matches:
            cansnps.append(CanSnp(column=column, bases=bases, missing=missing))
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
                missing=missing,
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
            f"genome with a base at the site, so which side carries the derived base cannot be "
            f"told without a genome outside both; each is reported against the first of them. "
            f"An outgroup genome in the alignment, placed at no node and called at the site, "
            f"settles it")
    if undecided:
        logging.warning(
            f"{undecided} SNPs fit both a node and a node below it, because the genomes that "
            f"tell the two apart have no base at the site; they are not reported as canSNPs")
    logging.info(f"Found {canonical} canSNPs")

    return cansnps


def ranking(score):
    """The order SNPs of a node are taken in: those in sequence that occurs once
    in the reference first, then those with fewer genomes missing a base, then
    the best total, then by contig and position, so that the same scores
    always give the same choice."""
    return (getattr(score, "copies", 1) > 1, getattr(score, "missing", 0), -score.total,
            score.contig, score.position)


def select_per_node(scores, top=None, spacing=0):
    """The best ``top`` SNPs of each node, at least ``spacing`` bases apart.

    The SNPs of a node are taken in order of :func:`ranking`, passing over any
    that lie within ``spacing`` bases of one already taken, so that those
    chosen for a node do not share a locus that one deletion or recombination
    event could remove. With ``top`` None every SNP so spaced is taken.

    Returns:
        For each node, the chosen SNPs, best first.
    """
    by_node = defaultdict(list)
    for score in scores:
        by_node[score.node].append(score)
    chosen = {}
    for node, node_scores in by_node.items():
        picked = []
        for score in sorted(node_scores, key=ranking):
            if top is not None and len(picked) >= top:
                break
            if any(other.contig == score.contig
                   and abs(other.position - score.position) < spacing for other in picked):
                continue
            picked.append(score)
        chosen[node] = picked
    return chosen


def apply_filters(scores, cutoff, top=None, spacing=0):
    """
    Mark the scores that fail the cutoff, or are not among the best of their node.

    Nothing is dropped: a score at or below the cutoff gains the FILTER code
    LowScore, and one that passes it but is not among the ``top`` best of its
    node gains NotTop, so that the output still shows what was considered. The
    best are taken in order of score, passing over any within ``spacing`` bases
    of one already taken for the node, so that the SNPs chosen for a node do not
    share a locus that a single deletion or recombination event could remove.

    Args:
        scores: Scored SNPs, each with a total, a node and a list of filters.
        cutoff: A score must be above this to pass.
        top: How many of the best scores of each node pass. All of them if None.
        spacing: The fewest bases between two SNPs taken for the same node.

    Returns:
        The same scores, in the order given.
    """
    for score in scores:
        if score.total <= cutoff and "LowScore" not in score.filters:
            score.filters.append("LowScore")

    if top is not None:
        candidates = [score for score in scores if not score.filters]  # not LowScore, Masked
        chosen = {id(score) for picked in select_per_node(candidates, top, spacing).values()
                  for score in picked}
        for score in candidates:
            if id(score) not in chosen:
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
    node_counts = dict.fromkeys(nodes_data, 0)

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
