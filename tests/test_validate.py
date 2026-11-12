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
    assert cansnp_trf.run_trf("in.fasta", str(tmp_path / "out.trf")) == 1
    assert "is not on the PATH" in caplog.text


def test_a_failing_trf_run_is_reported(monkeypatch, tmp_path, caplog):
    # A non-zero exit used to be ignored, leaving an empty report that the
    # evaluator then read as a genome with no repeats anywhere.
    class Completed:
        returncode = 1
        stderr = "trf: cannot read the input"

    monkeypatch.setattr(cansnp_trf.shutil, "which", lambda name: "/usr/bin/trf")
    monkeypatch.setattr(cansnp_trf.subprocess, "run", lambda *a, **k: Completed())
    assert cansnp_trf.run_trf("in.fasta", str(tmp_path / "out.trf")) == 1
    assert "must not be used" in caplog.text


def test_a_successful_trf_run_reports_success(monkeypatch, tmp_path):
    class Completed:
        returncode = 0
        stderr = ""

    monkeypatch.setattr(cansnp_trf.shutil, "which", lambda name: "/usr/bin/trf")
    monkeypatch.setattr(cansnp_trf.subprocess, "run", lambda *a, **k: Completed())
    assert cansnp_trf.run_trf("in.fasta", str(tmp_path / "out.trf")) == 0
