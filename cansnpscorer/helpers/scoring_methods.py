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
    whether the SNP falls in a gene or CDS, which makes its context better
    defined.
``mutation``
    whether the substitution is synonymous in the CDS it falls in, read in
    that CDS's frame, which makes it more likely to be neutral; one that
    introduces or removes a stop codon, or lies outside any CDS, scores 0.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from cansnpscorer.helpers.data_utils import codon_comparer, gc_content, neighbour_distance, search
from cansnpscorer.model import Contig
from cansnpscorer.model import complement as complement_base

BASES = frozenset("ACGT")

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
    """Bases covered by a gene or CDS on this contig as (start, end), sorted and joined."""
    cds: Sequence[Sequence[int]] = field(default_factory=list)
    """CDS intervals on this contig as (start, end, strand, phase), sorted by start."""
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


def gene_score(context: SnpContext, limit: int) -> float:
    """Whether the SNP falls in a gene or CDS.

    The score is 1 inside one and falls to 0 at the limit, so a SNP in a well
    defined coding context scores highest.
    """
    distance, _interval = search(context.position, context.genes, limit)
    return 1.0 - distance


def codons_at(context: SnpContext) -> Optional[List[str]]:
    """The codon each allele gives, read in the frame of the CDS the SNP is in.

    A CDS is read from its 5' end, skipping ``phase`` bases first, as GFF3
    defines the phase: forward from its start on the plus strand, and backward
    from its end on the minus strand, complemented. Returns None outside every
    CDS, for a SNP in the skipped bases, and for a codon the CDS ends inside.
    When CDS overlap, the first one that holds the SNP is read.
    """
    if len(context.alleles) != 2:
        return None
    sequence = context.contig.sequence
    position = context.position

    for start, end, strand, phase in context.cds:
        if start > position:
            break
        if position > end:
            continue
        if strand == "-":
            # Transcribed from the end: the first codon begins at end - phase
            # and each codon reads towards lower coordinates.
            first = end - phase
            if position > first:
                continue
            index = (first - position) % 3
            top = position + index
            positions = [top, top - 1, top - 2]
            if positions[-1] < start:
                continue
            codon = [complement_base(sequence[at]) for at in positions]
            alleles = [complement_base(allele) for allele in context.alleles]
        else:
            first = start + phase
            if position < first:
                continue
            index = (position - first) % 3
            positions = [position - index + offset for offset in range(3)]
            if positions[-1] > end:
                continue
            codon = [sequence[at] for at in positions]
            alleles = list(context.alleles)
        if not set(codon) <= BASES:
            return None
        codons = []
        for allele in alleles:
            changed = list(codon)
            changed[index] = allele
            codons.append("".join(changed))
        return codons
    return None


def mutation_score(context: SnpContext) -> float:
    """1 for a synonymous substitution, 0.5 for one that changes the amino acid,
    and 0 for one that gains or loses a stop codon or lies outside every CDS."""
    codons = codons_at(context)
    if codons is None:
        return 0.0
    return codon_comparer(codons, "+")


#: The INFO key each term is written under.
INFO_KEYS = {
    "neighbour_score": "NEIGHBOUR",
    "gc_score": "GC",
    "tandem_score": "TANDEM",
    "gene_score": "GENE",
    "mutation_score": "MUTATION",
}


@dataclass
class ScoredSnp:
    """One scored canSNP, addressed by contig and 1-based position."""

    contig: str
    position: int
    node: str
    total: int
    terms: Dict[str, int]
    derived: Optional[str] = None
    ancestral: Optional[str] = None
    filters: List[str] = field(default_factory=list)
    """The FILTER codes it failed: LowScore, NotTop or Masked. Empty when it passed."""

    @property
    def passed(self) -> bool:
        return not self.filters


def score_snp(context: SnpContext, tandem_limit: int = 50, gene_limit: int = 0,
              gc_goal: float = 0.5, weights: Optional[Dict[str, float]] = None) -> ScoredSnp:
    """Score one SNP.

    Every term is a percentage rounded to a whole number, and the total is
    their weighted sum.
    """
    weights = DEFAULT_WEIGHTS if weights is None else weights

    raw = {
        "neighbour_score": neighbour_score(context),
        "gc_score": gc_score(context, gc_goal),
        "tandem_score": repeat_score(context, tandem_limit),
        "gene_score": gene_score(context, gene_limit),
        "mutation_score": mutation_score(context),
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
                 repeats: Dict[str, Sequence[Sequence[int]]],
                 cds: Optional[Dict[str, Sequence[Sequence[int]]]] = None) -> List[SnpContext]:
    """Build the contexts for a set of SNPs, filling in each SNP's neighbours.

    Args:
        snps: (contig, 0-based position within the contig, node, alleles) per SNP.
        genes: Gene and CDS intervals per contig name.
        repeats: Tandem repeat intervals per contig name.
        cds: CDS intervals with their phase per contig name.

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
                cds=sorted((cds or {}).get(contig.name, []), key=lambda interval: interval[0]),
            )
        )
    return contexts
