"""Sampling SNP sites from several alignments into one matrix for a tree."""

import pytest
from conftest import cli_run

from cansnpscorer.model import Reference
from cansnpscorer.sites import sample_sites

REFERENCE = Reference.from_records([("c1", "ACGTACGT")])


def vcf_with_positions(tmp_path, *positions):
    path = tmp_path / "snps.vcf"
    lines = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    lines += [f"c1\t{position}\t.\t{REFERENCE.base(position - 1)}\tN\t.\tPASS\t." for position in positions]
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def two_alignments(tmp_path):
    """Alignment A without gaps; alignment B with a gap column after reference
    base 2, so that reference position 5 is column 5 of A and column 6 of B."""
    a_path = tmp_path / "a.fasta"
    a_path.write_text(">ref\nACGTACGT\n>genomeA1\nACGTTCGT\n>shared\nACGTCCGT\n")
    b_path = tmp_path / "b.fasta"
    b_path.write_text(">ref\nAC-GTACGT\n>genomeB1\nACAGTGCGT\n>shared\nACAGTACGT\n")
    return str(a_path), str(b_path)


def test_two_alignments_with_different_gaps_give_one_column_set(tmp_path):
    a_path, b_path = two_alignments(tmp_path)
    rows = sample_sites(vcf_with_positions(tmp_path, 5), REFERENCE, [a_path, b_path], sites=1)
    assert rows["genomeA1"] == "T" and rows["genomeB1"] == "G"
    assert rows["ref"] == "A"
    assert len(set(map(len, rows.values()))) == 1


def test_a_record_in_two_alignments_is_taken_from_the_first(tmp_path):
    a_path, b_path = two_alignments(tmp_path)
    rows = sample_sites(vcf_with_positions(tmp_path, 5), REFERENCE, [a_path, b_path], sites=1)
    assert rows["shared"] == "C"
    assert list(rows) == ["ref", "genomeA1", "shared", "genomeB1"]


def test_sampling_is_reproducible(tmp_path):
    a_path, b_path = two_alignments(tmp_path)
    vcf = vcf_with_positions(tmp_path, 1, 2, 4, 5, 7, 8)
    first = sample_sites(vcf, REFERENCE, [a_path, b_path], sites=3, seed=7)
    assert first == sample_sites(vcf, REFERENCE, [a_path, b_path], sites=3, seed=7)
    assert len(first["genomeA1"]) == 3


def test_the_sample_is_that_of_python_random_on_the_positions_in_file_order(tmp_path):
    # The benchmark's alltree.py took random.seed(0); random.sample(positions)
    # with the positions in VCF order, then sorted them.
    import random

    a_path, _b_path = two_alignments(tmp_path)
    positions = [8, 2, 5, 1, 7, 4]
    rows = sample_sites(vcf_with_positions(tmp_path, *positions), REFERENCE, [a_path], sites=3, seed=0)
    taken = sorted(random.Random(0).sample([p - 1 for p in positions], 3))
    assert rows["genomeA1"] == "".join("ACGTTCGT"[p] for p in taken)


def test_more_sites_than_snps_takes_every_snp(tmp_path):
    a_path, _b_path = two_alignments(tmp_path)
    rows = sample_sites(vcf_with_positions(tmp_path, 2, 5), REFERENCE, [a_path], sites=100)
    assert rows["genomeA1"] == "CT"


def test_a_ragged_alignment_is_refused(tmp_path):
    path = tmp_path / "ragged.fasta"
    path.write_text(">ref\nACGTACGT\n>short\nACGT\n")
    with pytest.raises(ValueError, match="not aligned"):
        sample_sites(vcf_with_positions(tmp_path, 5), REFERENCE, [str(path)], sites=1)


def test_sites_from_the_command_line(tmp_path):
    a_path, b_path = two_alignments(tmp_path)
    reference = tmp_path / "reference.fasta"
    reference.write_text(">c1\nACGTACGT\n")
    output = tmp_path / "sites.fasta"
    assert (
        cli_run(
            "sites",
            "--snps",
            vcf_with_positions(tmp_path, 5),
            "--reference",
            str(reference),
            "--alignment",
            a_path,
            "--alignment",
            b_path,
            "--sites",
            "1",
            "--seed",
            "0",
            "-o",
            str(output),
        )
        == 0
    )
    assert output.read_text() == ">ref\nA\n>genomeA1\nT\n>shared\nC\n>genomeB1\nG\n"


def test_sites_reports_a_contig_the_reference_lacks(tmp_path):
    a_path, _b_path = two_alignments(tmp_path)
    reference = tmp_path / "reference.fasta"
    reference.write_text(">c1\nACGTACGT\n")
    vcf = tmp_path / "other.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nc9\t5\t.\tA\tT\t.\tPASS\t.\n")
    assert (
        cli_run(
            "sites",
            "--snps",
            str(vcf),
            "--reference",
            str(reference),
            "--alignment",
            a_path,
            "-o",
            str(tmp_path / "out.fasta"),
        )
        == 1
    )


def test_a_name_repeated_within_one_alignment_is_refused(tmp_path):
    path = tmp_path / "dup.fasta"
    path.write_text(">ref\nACGTACGT\n>g1\nACGTTCGT\n>g1\nACGTCCGT\n")
    with pytest.raises(ValueError, match="'g1' is used more than once"):
        sample_sites(vcf_with_positions(tmp_path, 5), REFERENCE, [str(path)], sites=1)
