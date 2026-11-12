"""Scoring of a canSNP for its suitability as an assay target.

Every score is a function of one ``SnpContext``, which carries the contig the
SNP sits on, its position within that contig, its alleles, its neighbouring
SNPs and the gene and repeat intervals of that contig. Nothing here works in
concatenated genome coordinates, so no offset bookkeeping is needed and a SNP
is scored the same way whichever contig it lies on.

The five terms are, each on a scale of 0 to 1 before weighting:

``neighbour``
    distance to the closest neighbouring SNP, since a nearby SNP interferes
    with an assay.
``gc``
    how close the surrounding GC content is to the goal.
``repeat``
    distance to the closest tandem repeat.
``gene``
    whether the SNP falls in a gene, which makes its context better defined.
``mutation``
    whether the substitution is synonymous, which makes it more likely to be
    neutral; a substitution that introduces or removes a stop codon scores 0.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from cansnpscorer.helpers.data_utils import codon_comparer, gc_content, neighbour_distance, search
from cansnpscorer.model import Contig, reverse_complement

#: The score terms, in the order they are reported.
TERMS = ("neighbour_score", "gc_score", "tandem_score", "gene_score", "mutation_score")

#: Weights applied before the terms are summed. The gene term counts twice,
#: which is the behaviour the tool has had since the scoring was introduced.
DEFAULT_WEIGHTS: Dict[str, float] = {"gene_score": 2.0}

#: The highest total a SNP can reach with the default weights.
MAXIMUM_SCORE = round(100 * sum(DEFAULT_WEIGHTS.get(term, 1.0) for term in TERMS))


@dataclass
class SnpContext:
    """One SNP and everything needed to score it, in contig coordinates."""

    contig: Contig
    position: int
    """0-based position within the contig."""
    node: str
    alleles: Sequence[str]
    previous: Optional[int] = None
    """0-based position of the closest SNP to the left on this contig."""
    following: Optional[int] = None
    """0-based position of the closest SNP to the right on this contig."""
    genes: Sequence[Sequence[int]] = field(default_factory=list)
    """Gene intervals on this contig as (start, end, strand), sorted by start."""
    repeats: Sequence[Sequence[int]] = field(default_factory=list)
    """Tandem repeat intervals on this contig as (start, end), sorted by start."""

    @property
    def derived(self) -> Optional[str]:
        """The base carried below the node, when the polarity is known."""
        return self.alleles[0] if self.alleles else None

    @property
    def ancestral(self) -> Optional[str]:
        return self.alleles[1] if len(self.alleles) > 1 else None

    @property
    def global_position(self) -> int:
        return self.contig.start + self.position

    @property
    def reported_position(self) -> int:
        """The 1-based position within the contig, as schemes record it."""
        return self.position + 1


def neighbour_score(context: SnpContext, minimum_distance: int = 25) -> float:
    """Distance to the closest neighbouring SNP, normalized to minimum_distance."""
    left = context.previous if context.previous is not None else 0
    right = context.following if context.following is not None else context.contig.length
    return neighbour_distance(context.position, left, right, minimum_distance)


def gc_score(context: SnpContext, goal: float, distance: int = 50) -> float:
    """How close the GC content around the SNP is to the goal."""
    return gc_content(context.position, context.contig.sequence, goal, distance)


def repeat_score(context: SnpContext, limit: int) -> float:
    """Distance to the closest tandem repeat, normalized to limit."""
    score, _interval = search(context.position, context.repeats, limit)
    return score


def gene_score(context: SnpContext, limit: int) -> Tuple[float, Optional[Sequence[int]]]:
    """Whether the SNP falls in a gene, and the gene it falls in.

    The score is 1 inside a gene and falls to 0 at the limit, so a SNP in a
    well defined coding context scores highest. The interval is returned only
    when the SNP is inside the gene, where its reading frame applies.
    """
    distance, interval = search(context.position, context.genes, limit)
    return 1.0 - distance, interval


def mutation_score(context: SnpContext, interval: Optional[Sequence[int]]) -> float:
    """Compare the codons the alleles produce, for a SNP inside a gene.

    Returns 0 for a SNP outside a gene, and for one whose codon is incomplete
    because the annotated gene is not a whole number of codons.
    """
    if interval is None:
        return 0.0
    if len(context.alleles) != 2:
        return 0.0

    start, end, strand = interval[0], interval[1], interval[2]
    gene = context.contig.sequence[start:end + 1]
    offset = context.position - start
    if not 0 <= offset < len(gene):
        return 0.0

    alleles = list(context.alleles)
    if strand == "-":
        gene = reverse_complement(gene)
        offset = len(gene) - 1 - offset
        alleles = [reverse_complement(allele) for allele in alleles]

    codons = []
    codon_start = offset - offset % 3
    for allele in alleles:
        codon = list(gene[codon_start:codon_start + 3])
        if len(codon) != 3:
            return 0.0
        codon[offset - codon_start] = allele
        codons.append("".join(codon))

    # The strand is already accounted for, so compare the codons as they read.
    return codon_comparer(codons, "+")


@dataclass
class ScoredSnp:
    """One scored SNP, addressed by contig and 1-based position."""

    HEADER = (
        "Contig", "Position", "Group", "Derived", "Ancestral", "Score",
        "NeighborDistanceScore", "GCContentScore", "TandemRepeatProximityScore",
        "GeneProximityScore", "MutationTypeScore",
    )

    contig: str
    position: int
    node: str
    total: int
    terms: Dict[str, int]
    derived: Optional[str] = None
    ancestral: Optional[str] = None

    def row(self) -> List:
        """The record as it is written to the score file, following HEADER."""
        return [self.contig, self.position, self.node, self.derived or "", self.ancestral or "",
                self.total, *(self.terms[term] for term in TERMS)]


def score_snp(context: SnpContext, tandem_limit: int = 50, gene_limit: int = 0,
              gc_goal: float = 0.5, weights: Optional[Dict[str, float]] = None) -> ScoredSnp:
    """Score one SNP.

    Every term is a percentage rounded to a whole number, and the total is
    their weighted sum.
    """
    weights = DEFAULT_WEIGHTS if weights is None else weights

    gene, interval = gene_score(context, gene_limit)
    raw = {
        "neighbour_score": neighbour_score(context),
        "gc_score": gc_score(context, gc_goal),
        "tandem_score": repeat_score(context, tandem_limit),
        "gene_score": gene,
        "mutation_score": mutation_score(context, interval),
    }

    terms = {term: round(raw[term] * weights.get(term, 1.0) * 100) for term in TERMS}
    return ScoredSnp(
        contig=context.contig.name,
        position=context.reported_position,
        node=context.node,
        total=sum(terms.values()),
        terms=terms,
        derived=context.derived,
        ancestral=context.ancestral,
    )


def contexts_for(snps: Sequence[Tuple[Contig, int, str, Sequence[str]]],
                 genes: Dict[str, Sequence[Sequence[int]]],
                 repeats: Dict[str, Sequence[Sequence[int]]]) -> List[SnpContext]:
    """Build the contexts for a set of SNPs, filling in each SNP's neighbours.

    Args:
        snps: (contig, 0-based position within the contig, node, alleles) per SNP.
        genes: Gene intervals per contig name.
        repeats: Tandem repeat intervals per contig name.

    Returns:
        One context per SNP, ordered by contig and position. Neighbours are the
        adjacent SNPs on the same contig, whether or not they are canonical.
    """
    ordered = sorted(snps, key=lambda snp: (snp[0].start, snp[1]))

    positions: Dict[str, List[int]] = {}
    for contig, position, _node, _alleles in ordered:
        positions.setdefault(contig.name, []).append(position)

    contexts: List[SnpContext] = []
    seen: Dict[str, int] = {}
    for contig, position, node, alleles in ordered:
        index = seen.get(contig.name, 0)
        seen[contig.name] = index + 1
        on_contig = positions[contig.name]
        contexts.append(
            SnpContext(
                contig=contig,
                position=position,
                node=node,
                alleles=tuple(alleles),
                previous=on_contig[index - 1] if index > 0 else None,
                following=on_contig[index + 1] if index + 1 < len(on_contig) else None,
                genes=sorted(genes.get(contig.name, []), key=lambda interval: interval[0]),
                repeats=sorted(repeats.get(contig.name, []), key=lambda interval: interval[0]),
            )
        )
    return contexts
