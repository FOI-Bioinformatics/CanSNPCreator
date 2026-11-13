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

from collections.abc import Sequence
from dataclasses import dataclass, field

from cansnpscorer.helpers.data_utils import codon_comparer, gc_content, neighbour_distance, search
from cansnpscorer.model import Contig
from cansnpscorer.model import complement as complement_base

BASES = frozenset("ACGT")

#: The score terms, in the order they are reported.
TERMS = ("neighbour_score", "gc_score", "tandem_score", "gene_score", "mutation_score")

#: Weights applied before the terms are summed. The gene term counts twice,
#: which is the behaviour the tool has had since the scoring was introduced.
DEFAULT_WEIGHTS: dict[str, float] = {"gene_score": 2.0}

#: Weights for a scheme typed from genome alignments rather than by an assay.
#: The GC content around a SNP matters to a primer and not to an alignment, so
#: that term is left out, and the gene term counts once: a SNP in a gene is
#: more likely to be present in a new genome, which is worth one term, not two.
TYPING_WEIGHTS: dict[str, float] = {"gc_score": 0.0, "gene_score": 1.0}

#: The weights of each profile, and the cutoff it applies when none is given;
#: None applies no cutoff, so that every node keeps the canSNPs it has.
PROFILES: dict[str, tuple[dict[str, float], int | None]] = {
    "assay": (DEFAULT_WEIGHTS, 500),
    "typing": (TYPING_WEIGHTS, None),
}

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
    previous: int | None = None
    """0-based position of the closest SNP to the left on this contig."""
    following: int | None = None
    """0-based position of the closest SNP to the right on this contig."""
    genes: Sequence[Sequence[int]] = field(default_factory=list)
    """Bases covered by a gene or CDS on this contig as (start, end), sorted and joined."""
    cds: Sequence[Sequence[int]] = field(default_factory=list)
    """CDS intervals on this contig as (start, end, strand, phase), sorted by start."""
    repeats: Sequence[Sequence[int]] = field(default_factory=list)
    """Tandem repeat intervals on this contig as (start, end), sorted by start."""
    nearest: int | None = None
    """Distance to the nearest variation among the scheme's genomes, as extract
    measured it over every column where they differ; None when not known."""

    @property
    def derived(self) -> str | None:
        """The base carried below the node, when the polarity is known."""
        return self.alleles[0] if self.alleles else None

    @property
    def ancestral(self) -> str | None:
        return self.alleles[1] if len(self.alleles) > 1 else None

    @property
    def global_position(self) -> int:
        return self.contig.start + self.position

    @property
    def reported_position(self) -> int:
        """The 1-based position within the contig, as schemes record it."""
        return self.position + 1


def neighbour_score(context: SnpContext, flank: int = 50) -> float:
    """How clean the flanks of the SNP are, from 0 to 1 over ``flank`` bases.

    The distance is to the nearest variation among the scheme's genomes when
    extract measured it, which counts every column where they differ, and
    otherwise to the nearest SNP in the SNP file. A contig end counts as well,
    on purpose: an assay needs flanking sequence on both sides of its SNP,
    which a SNP near the end of a contig of a draft assembly does not have.
    """
    # Distances are differences of position, so a neighbouring site is 1 away,
    # and so is the first position past either end of the contig.
    if context.nearest is not None:
        to_ends = min(context.position + 1, context.contig.length - context.position)
        return min(context.nearest, to_ends, flank) / float(flank)
    left = context.previous if context.previous is not None else -1
    right = context.following if context.following is not None else context.contig.length
    return neighbour_distance(context.position, left, right, flank)


def gc_score(context: SnpContext, goal: float, distance: int = 50,
             tolerance: float = 0.1) -> float:
    """How close the GC content around the SNP is to the goal; full marks
    within ``tolerance`` of it, as a primer is designed in a range of GC."""
    return gc_content(context.position, context.contig.sequence, goal, distance, tolerance)


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


def codons_at(context: SnpContext) -> list[str] | None:
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
    terms: dict[str, int]
    derived: str | None = None
    ancestral: str | None = None
    filters: list[str] = field(default_factory=list)
    """The FILTER codes it failed: LowScore, NotTop or Masked. Empty when it passed."""
    missing: int = 0
    """Genomes of the alignment without a base at the SNP; fewer is ranked first."""
    copies: int = 1
    """Highest number of times a k-mer holding the SNP occurs in the reference, both
    strands counted; above 1 the SNP is ranked after those in unique sequence."""
    locus: str = ""
    """The gene the SNP lies in, as the first and last base of its interval on the
    contig; empty outside every gene. Two canSNPs of a node in one locus can be
    lost together, so a set is chosen from different loci where it can be."""
    discordant: int = 0
    """Genomes of a typing report in which the SNP disagreed with the other canSNPs
    of its node or with the genome's path; fewer is ranked first."""
    observed_genomes: int = 0
    """Genomes of the typing report the SNP was observed in; 0 when it was not
    in the scheme the report typed with. Where some canSNPs of a node were
    observed and found sound, they are ranked before those never observed."""
    observed_missing: int = 0
    """Genomes of a typing report without a base at the SNP; counted with ``missing``."""

    @property
    def passed(self) -> bool:
        return not self.filters


def score_snp(context: SnpContext, tandem_limit: int = 50, gene_limit: int = 0,
              gc_goal: float = 0.5, weights: dict[str, float] | None = None,
              flank: int = 50, gc_tolerance: float = 0.1) -> ScoredSnp:
    """Score one SNP.

    Every term is a percentage rounded to a whole number, and the total is
    their weighted sum.
    """
    weights = DEFAULT_WEIGHTS if weights is None else weights

    raw = {
        "neighbour_score": neighbour_score(context, flank),
        "gc_score": gc_score(context, gc_goal, tolerance=gc_tolerance),
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


def contexts_for(snps: Sequence[tuple[Contig, int, str, Sequence[str]]],
                 genes: dict[str, Sequence[Sequence[int]]],
                 repeats: dict[str, Sequence[Sequence[int]]],
                 cds: dict[str, Sequence[Sequence[int]]] | None = None) -> list[SnpContext]:
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

    positions: dict[str, list[int]] = {}
    for snp in ordered:
        positions.setdefault(snp[0].name, []).append(snp[1])

    # Sorted once per contig and shared by its SNPs, not copied for each, which
    # on a genome of a few thousand genes and some 10^5 SNPs is gigabytes.
    def by_start(intervals: dict[str, Sequence[Sequence[int]]]) -> dict[str, list]:
        return {name: sorted(found, key=lambda interval: interval[0])
                for name, found in intervals.items()}

    sorted_genes, sorted_repeats, sorted_cds = by_start(genes), by_start(repeats), by_start(cds or {})
    contexts: list[SnpContext] = []
    seen: dict[str, int] = {}
    for snp in ordered:
        contig, position, node, alleles = snp[:4]
        nearest = snp[4] if len(snp) > 4 else None
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
                genes=sorted_genes.get(contig.name, []),
                repeats=sorted_repeats.get(contig.name, []),
                cds=sorted_cds.get(contig.name, []),
                nearest=nearest,
            )
        )
    return contexts
