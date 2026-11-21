"""Audit of an existing CanSNPer2-data scheme.

The scheme audited here is the one the export subcommand writes from the
planted dataset, so a check that fires on it fires on a scheme this tool itself
produced. Each test then breaks one rule in the written files and asserts that
the audit names it.
"""

import os

import pytest
from conftest import REFERENCES, SNPS, TREE, cli_run, edit, info_of, vcf_lines, vcf_rewrite

from cansnpscorer import vcf
from cansnpscorer.cansnp_audit import PROBLEM, RECOMMENDATION, audit_scheme
from cansnpscorer.cansnp_audit import main as audit_main
from cansnpscorer.helpers.file_utils import reference_parser
from cansnpscorer.scheme import SchemeError, read_scheme


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
    assert audit_main(outdir, reference=dataset.reference) == 0
    assert audit_main(outdir, reference=dataset.reference, strict=True) == 2

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
    edit(outdir, SNPS, "\tTEST1\tc1\t201\t", "\tTEST2\tc1\t201\t")
    findings = findings_of(outdir, dataset)
    assert "unknown-genome" in categories(findings, PROBLEM)


# -- against the reference genome --------------------------------------------


def test_a_position_disagreeing_with_the_reference_is_reported(written, dataset):
    outdir, _scores = written
    edit(outdir, SNPS, "\t201\tC\tG\n", "\t201\tC\tT\n")
    findings = findings_of(outdir, dataset)
    reference_base = [finding for finding in findings if finding.category == "reference-base"]
    assert reference_base and reference_base[0].severity == PROBLEM
    assert reference_base[0].node == "A"


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
        base = reference.contig(fields[4]).base(int(fields[5]) - 1).upper()
        fields[6], fields[7] = sorted(set("ACGT") - {base})[:2]
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
    # than the one the scheme carries. REF must be the reference's base.
    header, records = vcf_lines(scores)
    current = next(line.split("\t") for line in records if info_of(line.split("\t"))
                   .get("NODE") == "A.1")
    ref = dataset.sequences["ref"][449]
    alt = next(base for base in "ACGT" if base != ref)
    info = ";".join(f"{key}=600" if key == "SCORE" else f"{key}={value}"
                    for key, value in info_of(current).items() if key != "AA")
    added = ["c1", "450", ".", ref, alt, ".", "PASS", f"{info};AA={ref}"]
    better = tmp_path / "better.vcf"
    better.write_text("\n".join(header + records + ["\t".join(added)]) + "\n")

    reference = reference_parser(dataset.reference)
    findings = findings_of(outdir, dataset, scores=vcf.read_scores(str(better), reference))
    offered = [finding for finding in findings if finding.category == "better-available"]
    assert [finding.node for finding in offered] == ["A.1"]
    assert "450" in offered[0].detail


def test_a_node_the_score_file_does_not_hold_is_reported(written, dataset, tmp_path):
    outdir, scores = written
    without = vcf_rewrite(scores, tmp_path / "without.vcf",
                          keep=lambda fields: info_of(fields).get("NODE") != "A.1")
    findings = findings_of(outdir, dataset, scores=vcf.read_scores(without))
    assert [finding.node for finding in findings if finding.category == "not-scored"] == ["A.1"]


def test_the_audit_runs_with_a_score_file(written, dataset):
    outdir, scores = written
    assert audit_main(outdir, reference=dataset.reference, scores=scores) == 0


def test_a_score_file_without_a_reference_is_refused(written, dataset):
    outdir, scores = written
    assert audit_main(outdir, scores=scores) == 1


# -- the command line --------------------------------------------------------


def test_the_audit_runs_from_the_command_line(written, dataset, tmp_path):
    outdir, _scores = written
    report = tmp_path / "findings.tsv"
    assert cli_run("audit", "--scheme-dir", outdir, "--organism", "testus_organismus",
                   "--reference", dataset.reference, "-o", str(report)) == 0
    assert "min-required-hits" in report.read_text()


def test_a_broken_scheme_exits_two_from_the_command_line(written, dataset):
    outdir, _scores = written
    edit(outdir, SNPS, "\tG\tC\n", "\t/\tC\n")
    assert cli_run("audit", "--scheme-dir", outdir,
                   "--reference", dataset.reference) == 2


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


def test_a_scheme_snp_in_a_masked_region_is_reported(written, dataset, tmp_path):
    outdir, _scores = written
    masks = tmp_path / "masks.bed"
    masks.write_text("c1\t395\t405\n")  # holds A.1's SNP at c1:401
    report = tmp_path / "findings.tsv"
    assert audit_main(outdir, reference=dataset.reference, masks=str(masks),
                      output=str(report)) == 0
    masked = [line for line in report.read_text().splitlines() if "\tmasked\t" in line]
    assert len(masked) == 1 and "\tA.1\t" in masked[0]
    assert audit_main(outdir, masks=str(masks)) == 1, "masks need the reference"
