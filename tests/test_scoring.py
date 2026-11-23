"""The five score terms, each named against the defect it guards.

The scorers used to take the SNP list, the concatenated genome, a correction
and a chromosome number, and recompute the contig boundary themselves at every
call. They now take one ``SnpContext`` in contig coordinates.
"""

import pytest

from cansnpscorer.helpers.scoring_methods import (
    SnpContext,
    contexts_for,
    gc_score,
    gene_score,
    mutation_score,
    neighbour_score,
    repeat_score,
    score_snp,
)
from cansnpscorer.model import Contig, Reference, reverse_complement

GENE = "ATGTTATAA"  # Met, Leu, stop


def context(sequence="ACGT" * 100, position=100, alleles=("A", "G"), **options):
    contig = options.pop("contig", None) or Contig(name="c1", sequence=sequence, start=0)
    return SnpContext(contig=contig, position=position, node="A.1", alleles=alleles, **options)


# -- the contig frame --------------------------------------------------------


def test_the_same_position_on_two_contigs_scores_the_same():
    # The wrappers indexed the reference list by a chromosome length rather
    # than a chromosome number, so anything beyond the first contig was scored
    # against the wrong sequence or raised IndexError.
    reference = Reference.from_records([("c1", "ACGT" * 100), ("c2", "ACGT" * 100)])
    genes = {"c1": [(50, 150, "+")], "c2": [(50, 150, "+")]}
    repeats = {"c1": [(10, 20)], "c2": [(10, 20)]}

    snps = [
        (reference.contig("c1"), 100, "A.1", ("A", "G")),
        (reference.contig("c2"), 100, "B.1", ("A", "G")),
    ]
    first, second = (score_snp(ctx) for ctx in contexts_for(snps, genes, repeats))
    assert first.terms == second.terms
    assert (first.contig, first.position) == ("c1", 101)
    assert (second.contig, second.position) == ("c2", 101)


def test_neighbours_are_taken_from_the_same_contig_only():
    reference = Reference.from_records([("c1", "ACGT" * 100), ("c2", "ACGT" * 100)])
    snps = [
        (reference.contig("c1"), 200, "A.1", ("A", "G")),
        (reference.contig("c2"), 210, "B.1", ("A", "G")),
    ]
    first, second = contexts_for(snps, {}, {})
    assert first.previous is None and first.following is None
    assert second.previous is None and second.following is None
    # Each is alone on its contig, so neither is held down by the other even
    # though they are ten bases apart in concatenated coordinates.
    assert neighbour_score(first) == pytest.approx(1.0)
    assert neighbour_score(second) == pytest.approx(1.0)


def test_a_snp_at_the_end_of_a_contig_has_no_flanking_sequence():
    reference = Reference.from_records([("c1", "ACGT" * 100)])
    at_end = contexts_for([(reference.contig("c1"), 399, "A.1", ("A", "G"))], {}, {})[0]
    assert neighbour_score(at_end) == pytest.approx(1 / 50)
    assert neighbour_score(at_end, flank=25) == pytest.approx(1 / 25)


def test_the_flank_is_judged_by_all_variation_when_extract_measured_it():
    reference = Reference.from_records([("c1", "ACGT" * 100)])
    snps = [(reference.contig("c1"), 200, "A.1", ("A", "G"), 10)]
    measured = contexts_for(snps, {}, {})[0]
    # Alone in the SNP file, but extract saw variation 10 bases away.
    assert measured.previous is None and measured.following is None
    assert neighbour_score(measured) == pytest.approx(10 / 50)
    # A contig end nearer than the variation still counts.
    near_end = contexts_for([(reference.contig("c1"), 2, "A.1", ("A", "G"), 40)], {}, {})[0]
    assert neighbour_score(near_end) == pytest.approx(3 / 50)
    far = contexts_for([(reference.contig("c1"), 200, "A.1", ("A", "G"), 1000)], {}, {})[0]
    assert neighbour_score(far) == pytest.approx(1.0)


def test_intervals_are_sorted_before_they_are_searched():
    # The scorers search the intervals with a binary search, which silently
    # misses on unsorted input; GFF and TRF files are not guaranteed to be
    # sorted.
    contig = Contig(name="c1", sequence="ACGT" * 100, start=0)
    unsorted_genes = {"c1": [(300, 350, "+"), (50, 150, "+")]}
    built = contexts_for([(contig, 100, "A.1", ("A", "G"))], unsorted_genes, {})[0]
    assert [interval[0] for interval in built.genes] == [50, 300]
    assert gene_score(built, 0) == 1.0


# -- the individual terms ----------------------------------------------------


def test_gc_content_near_a_contig_start_is_measured_on_a_clamped_window():
    # reference[snp - 50:snp + 50] wrapped on a negative start, which gave an
    # empty window and a ZeroDivisionError at the very beginning of a contig.
    at_start = context(sequence="GC" * 200, position=0)
    assert gc_score(at_start, 0.5, tolerance=0) == pytest.approx(0.0)
    balanced = context(sequence="AGCT" * 100, position=2)
    assert gc_score(balanced, 0.5) == pytest.approx(1.0)


def test_gc_content_within_the_tolerance_scores_as_the_goal():
    # 60% GC around the SNP: full marks within 0.1 of a goal of 0.5, and
    # beyond the band the term falls off as before from its edge.
    sixty = context(sequence="GGGCCCAAAT" * 40, position=200)
    assert gc_score(sixty, 0.5, tolerance=0.1) == pytest.approx(1.0)
    assert gc_score(sixty, 0.5, tolerance=0.0) == pytest.approx(0.8)
    assert gc_score(sixty, 0.5, tolerance=0.05) == pytest.approx(0.9)


def test_a_low_gc_goal_cannot_produce_a_negative_term():
    assert gc_score(context(sequence="GC" * 200), 0.2) == 0.0


def test_a_snp_near_a_gene_is_scored_without_crashing():
    # search returned a one element tuple when the SNP was outside the gene,
    # and the mutation scorer indexed element 1 of it, so any run with
    # --gene-limit greater than zero raised IndexError.
    near = context(position=90, sequence="ACGT" * 100)
    near.genes = [(100, 200, "+")]
    near.cds = [(100, 200, "+", 0)]
    score = gene_score(near, 200)
    assert 0.0 < score < 1.0, "close to a gene, but not in it"
    assert mutation_score(near) == 0.0, "no reading frame applies outside the CDS"
    assert score_snp(near, gene_limit=200).total >= 0


def test_a_snp_inside_a_gene_scores_full_proximity():
    inside = context(position=150)
    inside.genes = [(100, 200, "+")]
    assert gene_score(inside, 200) == 1.0


def test_repeat_proximity():
    with_repeat = context(position=100)
    with_repeat.repeats = [(80, 90)]
    assert repeat_score(with_repeat, 50) == pytest.approx(10 / 50)
    assert repeat_score(context(position=100), 50) == 1.0


# -- codons ------------------------------------------------------------------


def at(contig, position, alleles, cds):
    return SnpContext(contig=contig, position=position, node="A.1", alleles=alleles, cds=[cds])


def test_the_codon_is_read_in_the_frame_of_the_cds():
    contig = Contig(name="c1", sequence="TT" + GENE + "TT", start=0)
    cds = (2, 10, "+", 0)
    assert mutation_score(at(contig, 7, ("A", "G"), cds)) == 1.0, "TTA and TTG are both Leu"
    assert mutation_score(at(contig, 5, ("T", "G"), cds)) == 0.5, "Leu becomes Val"
    assert mutation_score(at(contig, 9, ("A", "T"), cds)) == 0.0, "the stop codon is lost"


def test_the_phase_moves_the_frame():
    # GFF3's phase is the number of bases to skip before the first codon. With
    # two leading bases in the CDS and a phase of 2, the codons are those of
    # GENE, as above.
    contig = Contig(name="c1", sequence="TT" + GENE + "TT", start=0)
    cds = (0, 10, "+", 2)
    assert mutation_score(at(contig, 7, ("A", "G"), cds)) == 1.0
    assert mutation_score(at(contig, 1, ("T", "G"), cds)) == 0.0, "before the first codon"


def test_a_cds_on_the_minus_strand_is_read_from_its_end():
    contig = Contig(name="c1", sequence=reverse_complement(GENE), start=0)
    cds = (0, 8, "-", 0)
    # The synonymous third base of codon two sits at offset 5 of the CDS,
    # which is position 3 of the contig on the minus strand.
    alleles = tuple(reverse_complement(base) for base in ("A", "G"))
    assert mutation_score(at(contig, 3, alleles, cds)) == 1.0
    # A phase of 1 on the minus strand skips the CDS's last base.
    shifted = Contig(name="c1", sequence=reverse_complement(GENE) + "A", start=0)
    assert mutation_score(at(shifted, 3, alleles, (0, 9, "-", 1))) == 1.0


def test_an_incomplete_codon_scores_zero():
    contig = Contig(name="c1", sequence=GENE, start=0)
    snp = at(contig, 6, ("T", "C"), (0, 7, "+", 0))  # the last codon is cut short
    assert mutation_score(snp) == 0.0


def test_a_gene_without_a_cds_gives_no_frame():
    # A gene feature says the SNP is in a gene; only a CDS says how it is read.
    snp = SnpContext(
        contig=Contig(name="c1", sequence=GENE, start=0),
        position=5,
        node="A.1",
        alleles=("T", "G"),
        genes=[(0, 8, "+")],
    )
    assert gene_score(snp, 0) == 1.0
    assert mutation_score(snp) == 0.0


def test_the_planted_dataset_scores_every_node(scored):
    nodes = {entry.node for entry in scored}
    assert nodes == {"root", "A", "B", "A.1", "A.2", "B.1", "B.2"}
    for entry in scored:
        assert entry.position >= 1, "positions are reported 1-based"
        assert entry.contig in {"c1", "c2"}
