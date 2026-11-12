"""Audit of an existing CanSNPer2-data scheme.

The scheme audited here is the one the export subcommand writes from the
planted dataset, so a check that fires on it fires on a scheme this tool itself
produced. Each test then breaks one rule in the written files and asserts that
the audit names it.
"""

import os
import sys

import pytest

from cansnpscorer.cansnp_audit import PROBLEM, RECOMMENDATION, audit_scheme
from cansnpscorer.cansnp_audit import main as audit_main
from cansnpscorer.cansnp_export import main as export_main
from cansnpscorer.cansnp_scorer import main as cli
from cansnpscorer.helpers.file_utils import read_score_file, reference_parser
from cansnpscorer.scheme import SchemeError, read_scheme

METADATA = dict(
    organism="testus_organismus",
    genome="TEST1",
    strain="testus",
    reference_label="Sjodin2026",
    genbank_id="GCA_000007125.1",
    refseq_id="GCF_000007125.1",
    assembly_name="ASM712v1",
)

TREE = "testus_organismus_tree.txt"
SNPS = "testus_organismus_snp.txt"
REFERENCES = "testus_organismus_references.txt"


def cli_run(*arguments):
    argv = sys.argv
    try:
        sys.argv = ["CanSNPCreator", *arguments]
        return cli()
    finally:
        sys.argv = argv


@pytest.fixture
def written(dataset, tmp_path):
    """A scheme written from the planted dataset, and its directory."""
    snp_file = str(tmp_path / "snps.tsv")
    scores = str(tmp_path / "scores.csv")
    outdir = str(tmp_path / "scheme")
    assert cli_run("extract", "--groupfile", dataset.groups_csv, "--hcfile", dataset.hc_yaml,
                   "--fastafile", dataset.alignment, "--reference-file", dataset.reference,
                   "-o", snp_file) == 0
    assert cli_run("evaluate", "--snp-file", snp_file, "--reference-file", dataset.reference,
                   "--gff-file", dataset.gff, "--trf-file", dataset.trf,
                   "--hcfile", dataset.hc_yaml, "--output-file", scores,
                   "--cutoff-value", "0") == 0
    assert export_main(scores, dataset.hc_yaml, dataset.reference, outdir, **METADATA) == 0
    return outdir, scores


def edit(directory, name, replace, with_text):
    path = os.path.join(directory, name)
    with open(path) as handle:
        text = handle.read()
    assert replace in text, replace
    with open(path, "w") as handle:
        handle.write(text.replace(replace, with_text))


def findings_of(directory, dataset=None, **arguments):
    reference = reference_parser(dataset.reference) if dataset else None
    return audit_scheme(read_scheme(directory), reference, **arguments)


def categories(findings, severity=None):
    return {finding.category for finding in findings
            if severity is None or finding.severity == severity}


# -- reading a scheme back in ------------------------------------------------


def test_the_written_scheme_is_read_back_whole(written):
    outdir, _scores = written
    scheme = read_scheme(outdir)
    assert scheme.organism == "testus_organismus"
    assert scheme.problems == []
    assert scheme.hierarchy is not None
    assert scheme.hierarchy.root == "root"
    assert len(scheme.snps) == len(scheme.nodes) == 7
    assert set(scheme.references_by_genome) == {"TEST1"}


def test_a_directory_with_no_scheme_is_refused(tmp_path):
    with pytest.raises(SchemeError, match="holds no"):
        read_scheme(str(tmp_path))


def test_two_schemes_in_one_directory_must_be_told_apart(written):
    outdir, _scores = written
    open(os.path.join(outdir, "other_tree.txt"), "w").close()
    with pytest.raises(SchemeError, match="--organism"):
        read_scheme(outdir)


# -- a scheme with nothing wrong ---------------------------------------------


def test_a_sound_scheme_reports_no_problem(written, dataset):
    outdir, _scores = written
    findings = findings_of(outdir, dataset)
    assert categories(findings, PROBLEM) == set()


def test_shallow_branches_are_reported_as_uncallable(written, dataset):
    # root, A and B have fewer than three SNPs on the path from the root, so
    # CanSNPer2 cannot call them with its default --min_required_hits.
    outdir, _scores = written
    findings = findings_of(outdir, dataset)
    hits = [finding for finding in findings if finding.category == "min-required-hits"]
    assert len(hits) == 1
    assert hits[0].severity == RECOMMENDATION
    assert hits[0].action == "extend"
    for node in ("root", "A", "B"):
        assert node in hits[0].node
    assert "A.1" not in hits[0].node


def test_a_lower_requirement_leaves_nothing_to_extend(written, dataset):
    outdir, _scores = written
    findings = findings_of(outdir, dataset, min_required_hits=1)
    assert "min-required-hits" not in categories(findings)


def test_the_exit_status_separates_problems_from_recommendations(written, dataset, tmp_path):
    outdir, _scores = written
    assert audit_main(outdir, reference_file=dataset.reference) == 0
    assert audit_main(outdir, reference_file=dataset.reference, strict=True) == 2

    report = tmp_path / "findings.tsv"
    assert audit_main(outdir, output=str(report)) == 0
    lines = report.read_text().splitlines()
    assert lines[0].split("\t") == ["severity", "action", "category", "node", "detail"]
    assert len(lines) > 1


# -- what the audit finds ----------------------------------------------------


def test_a_node_without_a_snp_is_reported(written, dataset):
    outdir, _scores = written
    edit(outdir, SNPS, "A.1\ttestus", "A.9\ttestus")
    findings = findings_of(outdir, dataset)
    assert {"no-snp", "unknown-node"} <= categories(findings, PROBLEM)


def test_a_node_with_two_snps_is_reported(written, dataset):
    outdir, _scores = written
    with open(os.path.join(outdir, SNPS)) as handle:
        lines = handle.read().splitlines()
    repeated = next(line for line in lines if line.startswith("A.1\t"))
    with open(os.path.join(outdir, SNPS), "a") as handle:
        handle.write(repeated.replace("\t401\t", "\t411\t") + "\n")

    findings = findings_of(outdir, dataset)
    assert "repeated-snp" in categories(findings, PROBLEM)
    assert "spacing" in categories(findings, RECOMMENDATION), "and they are 10 bases apart"


def test_a_node_under_two_parents_is_reported(written, dataset):
    outdir, _scores = written
    with open(os.path.join(outdir, TREE), "a") as handle:
        handle.write("root\tB\tA.1\n")
    findings = findings_of(outdir, dataset)
    assert "two-parents" in categories(findings, PROBLEM)


@pytest.mark.parametrize("replace,with_text,category", [
    # The three published Yersinia canSNPs carry a slash where the derived base
    # belongs, and can never be called derived.
    ("\tG\tC\n", "\t/\tC\n", "allele"),
    ("\tG\tC\n", "\tC\tC\n", "allele"),
    ("\t201\t", "\t401\t", "shared-position"),
    ("\t201\t", "\t0\t", "position"),
    ("\t201\t", "\t99999\t", "position"),
])
def test_an_unusable_snp_is_reported(written, dataset, replace, with_text, category):
    outdir, _scores = written
    edit(outdir, SNPS, replace, with_text)
    assert category in categories(findings_of(outdir, dataset), PROBLEM)


@pytest.mark.parametrize("replace,with_text,category", [
    ("GCA_000007125.1", "GCA_7125.1", "accession"),
    ("GCF_000007125.1", "GCF_7125.1", "accession"),
    ("TEST1\ttestus", "TES_T1\ttestus", "genome-name"),
])
def test_a_reference_row_that_breaks_the_rules_is_reported(written, dataset, replace, with_text,
                                                           category):
    outdir, _scores = written
    edit(outdir, REFERENCES, replace, with_text)
    assert category in categories(findings_of(outdir, dataset), PROBLEM)


def test_a_snp_on_an_unknown_genome_is_reported(written, dataset):
    outdir, _scores = written
    edit(outdir, SNPS, "\tTEST1\t201\t", "\tTEST2\t201\t")
    findings = findings_of(outdir, dataset)
    assert "unknown-genome" in categories(findings, PROBLEM)


# -- against the reference genome --------------------------------------------


def test_a_position_disagreeing_with_the_reference_is_reported(written, dataset):
    outdir, _scores = written
    edit(outdir, SNPS, "\t201\tC\tG\n", "\t201\tC\tT\n")
    findings = findings_of(outdir, dataset)
    reference_base = [finding for finding in findings if finding.category == "reference-base"]
    assert reference_base and reference_base[0].severity == PROBLEM
    assert "A" == reference_base[0].node


def test_coordinates_taken_against_another_assembly_are_reported(written, dataset):
    # A scheme whose positions were taken against another version of the
    # assembly presents itself as most of its SNPs sitting on a base that is
    # neither allele. Here every row is given two alleles that the reference
    # base at its position is not, which is that state exactly.
    outdir, _scores = written
    reference = reference_parser(dataset.reference)
    with open(os.path.join(outdir, SNPS)) as handle:
        header, *rows = handle.read().splitlines()

    changed = []
    for row in rows:
        fields = row.split("\t")
        base = reference.base(int(fields[4]) - 1).upper()
        fields[5], fields[6] = sorted(set("ACGT") - {base})[:2]
        changed.append("\t".join(fields))
    with open(os.path.join(outdir, SNPS), "w") as handle:
        handle.write("\n".join([header, *changed]) + "\n")

    findings = findings_of(outdir, dataset)
    assert "coordinates" in categories(findings, PROBLEM)
    assert sum(1 for finding in findings if finding.category == "reference-base") == len(rows)


def test_the_reference_base_is_not_checked_without_a_reference_file(written, dataset):
    outdir, _scores = written
    edit(outdir, SNPS, "\t201\tC\tG\n", "\t201\tC\tT\n")
    assert "reference-base" not in categories(findings_of(outdir))


# -- against scored candidates -----------------------------------------------


def test_a_better_scoring_candidate_is_offered(written, dataset, tmp_path):
    outdir, scores = written
    # A second candidate for A.1, elsewhere on the contig and scoring higher
    # than the one the scheme carries.
    with open(scores) as handle:
        rows = handle.read().splitlines()
    current = next(row for row in rows if row.split(",")[2] == "A.1")
    fields = current.split(",")
    fields[1], fields[5] = "450", "600"
    better = tmp_path / "better.csv"
    better.write_text("\n".join([*rows, ",".join(fields)]) + "\n")

    findings = findings_of(outdir, dataset, scores=read_score_file(str(better)))
    offered = [finding for finding in findings if finding.category == "better-available"]
    assert [finding.node for finding in offered] == ["A.1"]
    assert "450" in offered[0].detail


def test_a_node_the_score_file_does_not_hold_is_reported(written, dataset, tmp_path):
    outdir, scores = written
    with open(scores) as handle:
        rows = handle.read().splitlines()
    without = tmp_path / "without.csv"
    without.write_text("\n".join(row for row in rows if row.split(",")[2] != "A.1") + "\n")

    findings = findings_of(outdir, dataset, scores=read_score_file(str(without)))
    assert [finding.node for finding in findings if finding.category == "not-scored"] == ["A.1"]


def test_the_audit_runs_with_a_score_file(written, dataset):
    outdir, scores = written
    assert audit_main(outdir, reference_file=dataset.reference, score_file=scores) == 0


def test_a_score_file_without_a_reference_is_refused(written, dataset):
    outdir, scores = written
    assert audit_main(outdir, score_file=scores) == 1


# -- the command line --------------------------------------------------------


def test_the_audit_runs_from_the_command_line(written, dataset, tmp_path):
    outdir, _scores = written
    report = tmp_path / "findings.tsv"
    assert cli_run("audit", "--scheme-dir", outdir, "--organism", "testus_organismus",
                   "--reference-file", dataset.reference, "-o", str(report)) == 0
    assert "min-required-hits" in report.read_text()


def test_a_broken_scheme_exits_two_from_the_command_line(written, dataset):
    outdir, _scores = written
    edit(outdir, SNPS, "\tG\tC\n", "\t/\tC\n")
    assert cli_run("audit", "--scheme-dir", outdir,
                   "--reference-file", dataset.reference) == 2


# -- against a published scheme ----------------------------------------------


PUBLISHED = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "CanSNPer2-data", "source")


@pytest.mark.skipif(not os.path.isdir(PUBLISHED),
                    reason="CanSNPer2-data is not checked out beside this repository")
def test_the_published_yersinia_scheme_reports_its_known_alleles():
    # Three canSNPs of that scheme are published with a slash where the derived
    # base belongs, and are listed under known_issues in its schemes.toml. They
    # are what the allele check exists for, so the audit must find them without
    # being told.
    findings = audit_scheme(read_scheme(PUBLISHED, "yersinia_pestis"))
    alleles = [finding for finding in findings
               if finding.category == "allele" and finding.severity == PROBLEM]
    assert {finding.node for finding in alleles} == {"1.ANT2.b", "1.ORI2.c", "1.ORI3.b"}


@pytest.mark.skipif(not os.path.isdir(PUBLISHED),
                    reason="CanSNPer2-data is not checked out beside this repository")
@pytest.mark.parametrize("organism", ["francisella_tularensis", "bacillus_anthracis"])
def test_a_published_scheme_reads_without_a_problem(organism):
    scheme = read_scheme(PUBLISHED, organism)
    assert scheme.problems == []
    assert scheme.hierarchy is not None
    assert len(scheme.snps) == len(scheme.nodes)
    structural = [finding for finding in audit_scheme(scheme) if finding.severity == PROBLEM]
    assert structural == []
