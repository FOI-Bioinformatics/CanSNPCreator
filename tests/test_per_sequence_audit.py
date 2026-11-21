"""The audit reads each canSNP on its own reference sequence.

Since CanSNPer2 2.3.0 a scheme gives a canSNP as a sequence and a position on
it. The reference here holds two sequences of five bases, c1 of A and c2 of
C, so a canSNP read on the wrong sequence, or at a position counted along the
joined sequences, meets the wrong base.
"""

from cansnpscorer.cansnp_audit import (
    PROBLEM,
    check_against_reference,
    check_against_scores,
    check_alleles,
    check_masks,
    check_set_region,
    check_spacing,
)
from cansnpscorer.helpers.scoring_methods import ScoredSnp
from cansnpscorer.model import AlignmentMap, Reference
from cansnpscorer.scheme import Scheme, SchemeSnp

REFERENCE = Reference.from_records([("c1", "AAAAA"), ("c2", "CCCCC")])


def snp(node, sequence, position, derived="G", ancestral="C", rank=1):
    return SchemeSnp(node, "s", "r", "G1", sequence, position, derived, ancestral, rank=rank)


def scheme_with(*snps, snp_set=()):
    return Scheme(organism="t", snps=list(snps), snp_set=list(snp_set))


def categories(findings):
    return [finding.category for finding in findings]


def test_the_reference_base_is_read_on_the_snps_own_sequence():
    # c2:1 is C, the ancestral base; c1:1, where a joined count would read
    # position 1, is A.
    findings = check_against_reference(scheme_with(snp("B.2", "c2", 1)), REFERENCE)
    assert "reference-base" not in categories(findings)
    assert "position" not in categories(findings)


def test_a_reference_base_that_differs_names_the_sequence_and_position():
    findings = check_against_reference(scheme_with(snp("B.2", "c2", 3, "G", "T")), REFERENCE)
    (found,) = [finding for finding in findings if finding.category == "reference-base"]
    assert "'C' at c2:3" in found.detail
    assert "concatenated" not in found.detail


def test_a_position_past_the_end_of_its_sequence_is_reported():
    # Position 6 lies within the ten bases of the joined sequences, but not
    # within c1.
    findings = check_against_reference(scheme_with(snp("A", "c1", 6)), REFERENCE)
    (found,) = [finding for finding in findings if finding.category == "position"]
    assert found.severity == PROBLEM and "c1, which holds 5 bases" in found.detail


def test_a_sequence_the_reference_does_not_hold_is_reported():
    findings = check_against_reference(
        scheme_with(snp("A", "c1", 1, "A", "G"), snp("B", "CP000001.1", 2),
                    snp("C", "CP000001.1", 3)), REFERENCE)
    (found,) = [finding for finding in findings if finding.category == "sequence"]
    assert found.severity == PROBLEM
    assert "2 SNPs" in found.detail and "'CP000001.1'" in found.detail
    assert "c1, c2" in found.detail


def test_one_position_on_two_sequences_is_not_shared():
    assert check_alleles(scheme_with(snp("A", "c1", 3), snp("B", "c2", 3))) == []
    shared = check_alleles(scheme_with(snp("A", "c1", 3), snp("B", "c1", 3)))
    assert categories(shared) == ["shared-position"]
    assert "G1 c1:3" in shared[0].detail


def test_spacing_is_measured_along_one_sequence():
    assert check_spacing(scheme_with(snp("A", "c1", 3), snp("B", "c2", 4)), 50) == []
    (found,) = check_spacing(scheme_with(snp("A", "c2", 3), snp("B", "c2", 4)), 50)
    assert "on G1 c2" in found.detail


def test_a_set_is_clustered_only_along_one_sequence():
    spread = scheme_with(snp("A", "c1", 100), snp_set=[
        snp("A", "c1", 105, rank=2), snp("A", "c2", 102, rank=3)])
    assert check_set_region(spread, 10) == []
    together = scheme_with(snp("A", "c2", 100), snp_set=[
        snp("A", "c2", 105, rank=2), snp("A", "c2", 102, rank=3)])
    (found,) = check_set_region(together, 10)
    assert "G1 c2 (100-105)" in found.detail


def test_a_scored_candidate_is_matched_on_its_sequence():
    scheme = scheme_with(snp("B", "c2", 2))
    scores = [ScoredSnp(contig="c2", position=2, node="B", total=100, terms={}),
              ScoredSnp(contig="c2", position=4, node="B", total=200, terms={})]
    findings = check_against_scores(scheme, scores, REFERENCE)
    assert categories(findings) == ["better-available"]
    assert "scheme's SNP at c2:2 scores 100" in findings[0].detail
    assert "c2:4 scores 200" in findings[0].detail


def test_a_scored_candidate_on_another_sequence_is_not_the_schemes():
    scheme = scheme_with(snp("B", "c2", 2))
    scores = [ScoredSnp(contig="c1", position=2, node="B", total=100, terms={})]
    (found,) = check_against_scores(scheme, scores, REFERENCE)
    assert found.category == "not-scored" and "c2:2" in found.detail


def test_a_mask_holds_the_snps_on_its_own_sequence():
    masks = {"c2": [(0, 2)]}  # c2:1-3, 0-based and closed
    found = check_masks(scheme_with(snp("B", "c2", 2), snp("A", "c1", 2)), REFERENCE, masks)
    assert [finding.node for finding in found] == ["B"]
    assert "c2:2" in found[0].detail


def test_an_alignment_column_is_found_by_sequence_and_position():
    # The reference row has a gap in the first column.
    mapping = AlignmentMap("-AAAAACCCCC", REFERENCE)
    assert mapping.column_at("c1", 1) == 1
    assert mapping.column_at("c2", 1) == 6
    assert mapping.column_at("c2", 6) is None
    assert mapping.column_at("c3", 1) is None
