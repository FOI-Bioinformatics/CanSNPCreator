"""The checks that need more than one file to see, and the TRF runner."""

import pytest
from conftest import CLUSTERS

from cansnpscorer import cansnp_trf
from cansnpscorer.helpers.file_utils import (
    read_alignment,
    read_genome_map,
    read_hierarchy,
    reference_parser,
)
from cansnpscorer.validate import (
    check_alignment_against_reference,
    check_annotations,
    check_groups,
    report,
)


@pytest.fixture
def pieces(dataset):
    return (
        read_alignment(dataset.alignment),
        reference_parser(dataset.reference),
        read_hierarchy(dataset.hierarchy),
        read_genome_map(dataset.genome_map),
    )


# -- the alignment against the reference -------------------------------------


def test_consistent_inputs_report_nothing(pieces):
    alignment, reference, hierarchy, clusters = pieces
    assert check_alignment_against_reference(alignment, reference) == ([], [])
    problems, warnings = check_groups(alignment, clusters, hierarchy, "ref")
    assert (problems, warnings) == ([], [])


def test_a_reference_of_another_length_is_a_problem(pieces, tmp_path):
    alignment, reference, _hierarchy, _clusters = pieces
    shortened = tmp_path / "short.fasta"
    shortened.write_text(">c1\n" + reference.contigs[0].sequence[:100] + "\n")
    problems, _warnings = check_alignment_against_reference(
        alignment, reference_parser(str(shortened)))
    assert any("must be the reference genome" in problem for problem in problems)


def test_a_reference_of_the_same_length_but_other_sequence_is_a_problem(pieces, tmp_path):
    alignment, reference, _hierarchy, _clusters = pieces
    other = tmp_path / "other.fasta"
    other.write_text(">c1\n" + "A" * reference.length + "\n")
    problems, _warnings = check_alignment_against_reference(
        alignment, reference_parser(str(other)))
    assert any("differs from the reference genome" in problem for problem in problems)


def test_gaps_in_the_reference_row_are_a_warning(dataset, tmp_path):
    reference = reference_parser(dataset.reference)
    gapped = dict(dataset.sequences)
    gapped["ref"] = "-" + gapped["ref"]
    _problems, warnings = check_alignment_against_reference(gapped, reference)
    assert any("carries 1 gaps" in warning for warning in warnings)


# -- groups, hierarchy and alignment -----------------------------------------


def test_an_isolate_missing_from_the_alignment_is_a_problem(pieces):
    alignment, _reference, hierarchy, clusters = pieces
    clusters = {**clusters, "B.2": ["s4", "s99"]}
    problems, _warnings = check_groups(alignment, clusters, hierarchy, "ref")
    assert any("the alignment does not hold" in problem for problem in problems)


def test_a_cluster_missing_from_the_hierarchy_is_a_problem(pieces):
    alignment, _reference, hierarchy, clusters = pieces
    clusters = {**clusters, "C.1": ["s1"]}
    problems, _warnings = check_groups(alignment, clusters, hierarchy, "ref")
    assert any("the hierarchy does not hold" in problem for problem in problems)


def test_an_isolate_in_two_clusters_is_a_problem(pieces):
    alignment, _reference, hierarchy, clusters = pieces
    clusters = {**clusters, "B.1": ["s3", "s1"]}
    problems, _warnings = check_groups(alignment, clusters, hierarchy, "ref")
    assert any("in more than one cluster" in problem for problem in problems)


def test_a_genome_in_no_cluster_is_a_warning(pieces):
    # It still counts against every canSNP, which is why it is worth saying.
    alignment, _reference, hierarchy, _clusters = pieces
    clusters = {name: list(isolates) for name, isolates in CLUSTERS.items() if name != "B.2"}
    _problems, warnings = check_groups(alignment, clusters, hierarchy, "ref")
    assert any("are in no cluster" in warning for warning in warnings)
    assert any("have no isolates" in warning for warning in warnings)


def test_the_reference_genome_is_not_expected_to_be_grouped(pieces):
    alignment, _reference, hierarchy, clusters = pieces
    _problems, warnings = check_groups(alignment, clusters, hierarchy, "ref")
    assert not any("ref" in warning for warning in warnings)


# -- annotations -------------------------------------------------------------


def test_an_annotation_outside_its_contig_is_a_problem(dataset):
    reference = reference_parser(dataset.reference)
    problems, _warnings = check_annotations({"c2": [(10, 5000, "+")]}, reference, "gene")
    assert any("lie outside the contig" in problem for problem in problems)


def test_a_contig_without_annotations_is_a_warning(dataset):
    reference = reference_parser(dataset.reference)
    _problems, warnings = check_annotations({"c1": [(10, 20, "+")]}, reference, "gene")
    assert any("have no gene annotation" in warning for warning in warnings)


def test_report_says_whether_the_run_may_go_ahead():
    assert report([], ["a warning"]) is True
    assert report(["a problem"], []) is False


# -- the TRF runner ----------------------------------------------------------


def test_a_missing_trf_binary_is_reported(monkeypatch, tmp_path, caplog):
    monkeypatch.setattr(cansnp_trf.shutil, "which", lambda name: None)
    assert cansnp_trf.run_trf("in.fasta", str(tmp_path / "out.bed")) == 1
    assert "is not on the PATH" in caplog.text


NGS = """@c1 Testus organismus c1, complete sequence
60 95 4 9.0 4 90 0 50 25 25 25 25 1.9 ACGT ACGTACGTACGT
10 20 2 5.5 2 100 0 22 50 50 0 0 1.0 AC ACACACACACA
@c2 plasmid
120 160 8 5.1 8 85 5 40 25 25 25 25 2.0 ACGTACGT ACGTACGTACGT
"""


def fake_trf(monkeypatch, tmp_path, output=NGS, status=0):
    """A trf on the PATH that prints this output and exits with this status."""
    import stat
    import sys

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    script = bin_dir / "trf"
    script.write_text(f"#!{sys.executable}\nimport sys\n"
                      f"sys.stdout.write({output!r})\n"
                      f"sys.stderr.write('trf: cannot read the input' if {status} else '')\n"
                      f"sys.exit({status})\n")
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", f"{bin_dir}:" + __import__("os").environ["PATH"])


def test_trf_output_is_written_as_bed(monkeypatch, tmp_path):
    fake_trf(monkeypatch, tmp_path)
    output = tmp_path / "repeats.bed"
    assert cansnp_trf.run_trf("in.fasta", str(output)) == 0
    # 0-based and half open, sorted within each sequence, named by the ID alone.
    assert output.read_text().splitlines() == [
        "c1\t9\t20\tperiod=2,copies=5.5",
        "c1\t59\t95\tperiod=4,copies=9.0",
        "c2\t119\t160\tperiod=8,copies=5.1",
    ]


def test_a_failing_trf_run_writes_nothing(monkeypatch, tmp_path, caplog):
    # A non-zero exit used to be ignored, leaving an empty report that the
    # evaluator then read as a genome with no repeats anywhere.
    fake_trf(monkeypatch, tmp_path, output="", status=1)
    output = tmp_path / "repeats.bed"
    assert cansnp_trf.run_trf("in.fasta", str(output)) == 1
    assert "nothing was written" in caplog.text and "cannot read the input" in caplog.text
    assert not output.exists()


def test_a_repeat_before_any_sequence_is_refused():
    import pytest

    with pytest.raises(ValueError, match="before any sequence"):
        cansnp_trf.parse_ngs("60 95 4 9.0\n")


# -- regions in BED and GFF3 ---------------------------------------------------


def test_bed_and_gff3_regions_are_read_alike(tmp_path):
    from cansnpscorer.helpers.file_utils import read_intervals

    bed = tmp_path / "r.bed"
    bed.write_text("track name=repeats\n# comment\nc1\t9\t20\tx\nc1\t15\t30\ty\nc2\t0\t5\n")
    gff = tmp_path / "r.gff"
    gff.write_text("##gff-version 3\nc1\tGubbins\tCDS\t10\t20\t.\t+\t0\tID=a\n"
                   "c1\tGubbins\tCDS\t16\t30\t.\t+\t0\tID=b\nc2\tGubbins\tx\t1\t5\t.\t.\t.\t.\n")
    # Overlapping regions are joined, so a position is in a region or not.
    expected = {"c1": [(9, 29)], "c2": [(0, 4)]}
    assert read_intervals(str(bed)) == expected
    assert read_intervals(str(gff)) == expected


def test_a_malformed_region_is_refused(tmp_path):
    import pytest

    from cansnpscorer.helpers.file_utils import read_intervals

    bed = tmp_path / "r.bed"
    bed.write_text("c1\t20\t10\n")
    with pytest.raises(ValueError, match="impossible interval"):
        read_intervals(str(bed))
