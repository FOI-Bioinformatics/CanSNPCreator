"""evaluate counts how often the sequence around a canSNP occurs in the reference."""

import random

import pytest
from conftest import vcf_records
from test_cli import evaluate, extract

from cansnpscorer import vcf
from cansnpscorer.helpers.data_utils import ranking, select_per_node
from cansnpscorer.helpers.kmers import copy_numbers
from cansnpscorer.helpers.scoring_methods import ScoredSnp
from cansnpscorer.model import Reference

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def random_sequence(length, seed):
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(length))


REPEAT = random_sequence(60, 1)


def reference_with_repeat():
    """c1 holds the repeat at 100, and c2 holds its reverse complement at 200."""
    c1 = random_sequence(100, 2) + REPEAT + random_sequence(240, 3)
    c2 = random_sequence(200, 4) + REPEAT.translate(COMPLEMENT)[::-1] + random_sequence(140, 5)
    return Reference.from_records([("c1", c1), ("c2", c2)])


def test_a_site_in_unique_sequence_has_one_copy():
    assert copy_numbers(reference_with_repeat(), [("c1", 300)]) == {("c1", 300): 1}


def test_a_repeat_on_the_other_strand_of_another_contig_is_counted():
    reference = reference_with_repeat()
    copies = copy_numbers(reference, [("c1", 130), ("c2", 230)])
    assert copies == {("c1", 130): 2, ("c2", 230): 2}


def test_a_site_at_the_edge_of_a_repeat_is_counted_when_a_kmer_holding_it_lies_inside():
    reference = reference_with_repeat()
    # The last base of the repeat: the 31-mer ending there lies wholly inside it.
    assert copy_numbers(reference, [("c1", 159)])[("c1", 159)] == 2
    # One base past it, no 31-mer holding the site lies inside the repeat.
    assert copy_numbers(reference, [("c1", 160)])[("c1", 160)] == 1


def test_a_larger_step_still_finds_a_repeat_long_enough():
    # The repeat is 60 bases, at least k + step - 1 for k 31 and step 10.
    assert copy_numbers(reference_with_repeat(), [("c1", 130)], step=10)[("c1", 130)] == 2


def test_kmers_with_a_gap_base_are_not_counted():
    sequence = random_sequence(50, 6) + "N" * 40 + random_sequence(50, 7) + "N" * 40 + "ACGT"
    reference = Reference.from_records([("c1", sequence)])
    # Every 31-mer holding position 70 holds an N, and the N runs repeat.
    assert copy_numbers(reference, [("c1", 70)]) == {("c1", 70): 0}


def test_a_kmer_that_is_its_own_reverse_complement_is_counted_once_per_occurrence():
    # Only an even k allows it: ACGT reads the same on both strands. With
    # this step the one k-mer looked up for position 1 starts there.
    reference = Reference.from_records([("c1", "GACGTTG")])
    assert copy_numbers(reference, [("c1", 1)], k=4, step=100) == {("c1", 1): 1}
    # Two occurrences are two copies, not four.
    reference = Reference.from_records([("c1", "GACGTTGACGTG")])
    assert copy_numbers(reference, [("c1", 1)], k=4, step=100) == {("c1", 1): 2}


def test_size_and_step_must_be_positive():
    with pytest.raises(ValueError):
        copy_numbers(reference_with_repeat(), [("c1", 10)], k=0)


def scored(node, total, copies=1, position=100, missing=0):
    return ScoredSnp(
        contig="c1", position=position, node=node, total=total, terms={}, copies=copies, missing=missing
    )


def test_a_snp_in_repeated_sequence_is_ranked_after_one_in_unique_sequence():
    unique = scored("A", total=300, position=500)
    repeated = scored("A", total=600, copies=2, position=100)
    assert ranking(unique) < ranking(repeated)
    assert select_per_node([repeated, unique], top=1)["A"] == [unique]


def test_copies_come_before_missing_data_in_the_ranking():
    repeated_complete = scored("A", total=600, copies=3, position=100)
    unique_incomplete = scored("A", total=300, missing=2, position=500)
    assert ranking(unique_incomplete) < ranking(repeated_complete)


def test_evaluate_writes_the_copy_number_and_export_reads_it(monkeypatch, dataset, tmp_path):
    snps = str(tmp_path / "snps.vcf")
    assert extract(monkeypatch, dataset, snps) == 0
    scores = str(tmp_path / "scores.vcf")
    assert evaluate(monkeypatch, dataset, snps, scores) == 0
    records = vcf_records(scores)
    assert records and all(record.info["COPIES"] == "1" for record in records)
    read = vcf.read_scores(scores)
    assert all(score.copies == 1 for score in read)


def test_copy_kmer_zero_leaves_the_count_out(monkeypatch, dataset, tmp_path):
    snps = str(tmp_path / "snps.vcf")
    assert extract(monkeypatch, dataset, snps) == 0
    scores = str(tmp_path / "scores.vcf")
    assert evaluate(monkeypatch, dataset, snps, scores, "--copy-kmer", "0") == 0
    assert all("COPIES" not in record.info for record in vcf_records(scores))
    # A file without the field reads as unique sequence.
    assert all(score.copies == 1 for score in vcf.read_scores(scores))
