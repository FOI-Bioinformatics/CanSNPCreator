"""What the audit reports and how: the typing tables, the closing summary, the
findings that were corrected in the fourth audit round, and the Markdown
report."""

import os

from conftest import PLANTED, cli_run
from test_audit_round4 import chain_scheme, typed
from test_review_fixes import scheme_of
from test_typing_report import changed, run_report, table

from cansnpscorer.alignment_audit import check_against_alignment
from cansnpscorer.cansnp_audit import (
    check_against_reference,
    check_set_region,
    check_spacing,
    summarize,
)
from cansnpscorer.findings import MODIFY, problem, recommendation
from cansnpscorer.model import Reference
from cansnpscorer.report import markdown
from cansnpscorer.scheme import SchemeSnp
from cansnpscorer.typing_report import (
    GENOME_COLUMNS,
    NODE_COLUMNS,
    RANK_COLUMNS,
    SNP_COLUMNS,
    read_observed,
    write_report,
)

# -- the typing tables -------------------------------------------------------


def summary_of(report):
    return {row[0]: row[1:] for row in report.summary}


def test_the_reference_and_the_outgroups_are_typed_but_not_counted():
    alignment = {"ref": "AAAAAAAA", "g": "AGAGAGAA", "none": "AAAAAAAA", "out": "AAAAAAAA"}
    genomes, report = typed(chain_scheme(), alignment, outgroups=["out"])
    assert (genomes["ref"].role, genomes["out"].role, genomes["g"].role) == (
        "reference", "outgroup", "isolate")
    summary = summary_of(report)
    assert summary["records of the alignment"][0] == "4"
    assert summary["genomes"][0] == "2", "g and none"
    assert summary["typed"] == ("1", "50.0", "1", "50.0")
    assert summary["outgroups, not counted below"][0] == "1"
    assert summary["outgroups given a call"] == ("0", "NA", "0", "NA")


def test_the_summary_tallies_why_the_single_rule_stopped_and_names_the_nodes():
    # N at the canSNP of A in two genomes, ancestral there in a third whose
    # set of three still calls A derived.
    extra = [SchemeSnp("A", "", "", "G", 7, "G", "A", rank=2),
             SchemeSnp("A", "", "", "G", 8, "G", "A", rank=3)]
    alignment = {"ref": "AAAAAAAA", "m1": "AGANAGGG", "m2": "AGANAGGG", "a1": "AGAAAGGG"}
    scheme = chain_scheme(extra)
    scheme.snps = [row for row in scheme.snps if row.node != "B"]
    _genomes, report = typed(scheme, alignment)
    summary = summary_of(report)
    assert summary["single rule stopped at a canSNP without a base"][:2] == ("2", "66.7")
    assert summary["single rule stopped at an ancestral canSNP"][0] == "1"
    assert summary["nodes where the single rule stopped most often"][0] == (
        "A (2 missing); A (1 ancestral)")
    assert report.stops == [("A", "missing", 2), ("A", "ancestral", 1)]


def test_a_node_on_another_reference_is_named_as_the_reason(caplog):
    # A is annotated on genome H only, so in an alignment against G it has no
    # base, and the single rule stops there for a reason of the scheme's.
    scheme = scheme_of([["root"], ["root", "A"], ["root", "A", "A.1"]],
                       [("root", "G", 2, "G", "A"), ("A", "H", 4, "G", "A"),
                        ("A.1", "G", 6, "G", "A")])
    genomes, _report = typed(scheme, {"ref": "AAAAAAAA", "g": "AGAGAGAA"}, genome="G")
    assert genomes["g"].set_call == "A.1" and genomes["g"].single_call is None
    assert genomes["g"].single_stop == ("A", "other-reference")
    assert "other reference genomes only" in caplog.text


def test_nodes_no_genome_passes_through_are_listed():
    _genomes, report = typed(chain_scheme(), {"ref": "AAAAAAAA", "g": "AGAGAGAA"})
    assert report.unreached == ["B"]
    assert summary_of(report)["nodes no genome passes through"][2] == "1"


def test_the_ranks_table_sums_the_cansnps_of_each_rank():
    extra = [SchemeSnp("A", "", "", "G", 7, "G", "A", rank=2)]
    alignment = {"ref": "AAAAAAAA", "g": "AGAGAGNA", "h": "AGAGAGGA"}
    _genomes, report = typed(chain_scheme(extra), alignment)
    ranks = {row[0]: dict(zip(RANK_COLUMNS, row, strict=True)) for row in report.rank_rows}
    assert ranks["1"]["snps"] == "4" and ranks["2"]["snps"] == "1"
    assert (ranks["2"]["genomes"], ranks["2"]["missing"], ranks["2"]["missing_pct"]) == (
        "3", "1", "33.3")
    assert ranks["1"]["missing"] == "0"


def test_the_snp_table_gives_the_contig_and_no_field_is_left_empty(tmp_path):
    alignment = {"ref": "AAAAAAAA", "g": "AGAGAGAA", "none": "AAAAAAAA"}
    scheme = chain_scheme()
    reference = Reference.from_records([("c1", "AAAAA"), ("c2", "AAA")])
    from cansnpscorer.typing_report import type_alignment

    report = type_alignment(scheme, alignment, reference)
    prefix = str(tmp_path / "t")
    written = write_report(report, prefix, [("alignment", "a.fasta")])
    assert [os.path.basename(path) for path in written] == [
        "t_genomes.tsv", "t_snps.tsv", "t_nodes.tsv", "t_ranks.tsv", "t_summary.tsv",
        "t_run.tsv"]
    snps = {row["node"]: row for row in table(prefix + "_snps.tsv")}
    assert list(snps["A"]) == list(SNP_COLUMNS)
    assert (snps["A"]["contig"], snps["A"]["contig_position"]) == ("c1", "4")
    assert (snps["A.1"]["contig"], snps["A.1"]["contig_position"]) == ("c2", "1")
    for suffix, columns in (("_genomes.tsv", GENOME_COLUMNS), ("_nodes.tsv", NODE_COLUMNS),
                            ("_snps.tsv", SNP_COLUMNS)):
        for row in table(prefix + suffix):
            assert list(row) == list(columns)
            # derived_set is empty for a canSNP derived in no genome, and
            # set_apart for a node without a set of two canSNPs.
            assert all(value != "" for column, value in row.items()
                       if column not in ("derived_set", "set_apart")), (suffix, row)
    none = {row["genome"]: row for row in table(prefix + "_genomes.tsv")}["none"]
    assert (none["set_call"], none["set_depth"], none["single_stop_reason"]) == (
        "NA", "NA", "NA")
    assert table(prefix + "_run.tsv") == [{"item": "alignment", "value": "a.fasta"}]


def test_a_table_written_before_the_columns_were_renamed_is_still_read(tmp_path):
    old = ("node", "rank", "genome", "position", "derived_base", "ancestral_base", "genomes",
           "called", "missing", "derived", "ancestral", "other", "discordant_in_set",
           "ancestral_on_path", "derived_off_path")
    path = tmp_path / "old_snps.tsv"
    path.write_text("\t".join(old) + "\n"
                    + "\t".join(["A", "1", "G", "4", "G", "A", "9", "8", "1", "5", "3", "0",
                                 "1", "0", "2"]) + "\n")
    found = read_observed(str(path))[("A", 4)]
    assert (found.genomes, found.missing, found.discordant) == (8, 1, 2)


def test_many_genomes_and_cansnps_are_held_as_one_string_each():
    genomes, _report = typed(chain_scheme(), {"ref": "AAAAAAAA", "g": "AGAGAGAA"})
    assert genomes["g"].states == "DDDA"


def test_the_run_is_written_beside_the_tables(written, dataset, tmp_path):
    status, prefix = run_report(written, dataset, tmp_path)
    assert status in (0, 2)
    run = {row["item"]: row["value"] for row in table(prefix + "_run.tsv")}
    assert run["alignment"] == dataset.alignment and run["reference"] == dataset.reference
    assert run["min_required_hits"] == "3" and run["cansnpcreator"]


# -- the findings ------------------------------------------------------------


def test_spacing_compares_the_snps_cansnper2_calls_and_no_further_rank():
    scheme = scheme_of([["R"], ["R", "A"], ["R", "B"]],
                       [("A", "G", 100, "T", "G"), ("B", "G", 5000, "T", "G")])
    scheme.snp_set = [SchemeSnp("A", "", "", "G", 100, "T", "G", rank=1),
                      SchemeSnp("B", "", "", "G", 110, "T", "G", rank=2)]
    assert check_spacing(scheme, 50) == []
    scheme.snps.append(SchemeSnp("C", "s", "r", "G", 5010, "T", "G"))
    (found,) = check_spacing(scheme, 50)
    assert found.node == "B, C" and "assay of the other" in found.detail
    assert "primer" not in found.detail


def test_the_summary_counts_ranked_cansnps_and_names_each_node_once():
    scheme = scheme_of([["R"], ["R", "A"], ["R", "B"]],
                       [("A", "G", 100, "T", "G"), ("B", "G", 110, "T", "G")])
    scheme.snp_set = [SchemeSnp("A", "", "", "G", 100, "T", "G", rank=1),
                      SchemeSnp("A", "", "", "G", 900, "T", "G", rank=2)]
    findings = [recommendation(MODIFY, "spacing", "A, B", "close"),
                problem(MODIFY, "homoplasy", "A", "x"),
                problem(MODIFY, "homoplasy", "B", "y")]
    lines = summarize(scheme, findings)
    assert "2 SNPs, 1 further ranked canSNPs and 0 reference genomes" in lines[0]
    assert "Problems by category: homoplasy 2" in lines
    assert "Recommendations by category: spacing 1" in lines
    assert "Nodes with a SNP to modify (2): A, B" in lines


def wide_scheme_and_alignment(genomes=40):
    """root with the children A and B, twenty genomes below each, and an
    outgroup; A's SNP at position 2, B's at 4, A's second at 6."""
    scheme = scheme_of([["root"], ["root", "A"], ["root", "B"]],
                       [("A", "G", 2, "G", "A"), ("B", "G", 4, "G", "A")])
    alignment = {"ref": "AAAAAAAA", "out": "AAAAAAAA"}
    clusters = {"A": [], "B": []}
    for number in range(genomes):
        node = "A" if number % 2 == 0 else "B"
        name = f"{node.lower()}{number}"
        alignment[name] = "AGAAAGAA" if node == "A" else "AAAGAAAA"
        clusters[node].append(name)
    return scheme, alignment, clusters


def test_a_cansnp_missing_in_few_genomes_is_counted_and_not_listed_on_its_own():
    scheme, alignment, clusters = wide_scheme_and_alignment()
    alignment["a0"] = "ANAAAGAA"  # one of 42 genomes lacks A's SNP
    reference = Reference.from_records([("c1", alignment["ref"])])
    findings = check_against_alignment(scheme, alignment, reference, clusters)
    missing = [finding for finding in findings if finding.category == "missing-data"]
    assert [finding.node for finding in missing] == [""], "rolled up, not per canSNP"
    assert "1 further canSNPs" in missing[0].detail and "A (1)" in missing[0].detail
    (by_genome,) = [f for f in findings if f.category == "missing-data-genomes"]
    assert "a0 (1)" in by_genome.detail

    for name in clusters["A"][:5]:  # five of 42 is more than a twentieth
        alignment[name] = "ANAAAGAA"
    findings = check_against_alignment(scheme, alignment, reference, clusters)
    missing = [finding for finding in findings if finding.category == "missing-data"]
    assert [finding.node for finding in missing] == ["A"], "often missing: reported itself"


def test_single_site_says_whether_the_site_is_the_schemes_and_counts_its_set():
    scheme, alignment, clusters = wide_scheme_and_alignment(8)
    reference = Reference.from_records([("c1", alignment["ref"])])
    findings = check_against_alignment(scheme, alignment, reference, clusters)
    single = {f.node: f.detail for f in findings if f.category == "single-site"}
    assert "one of the scheme's own canSNPs" in single["B"]
    assert "nothing else in these genomes supports it" in single["B"]
    assert "A" not in single, "two sites separate A's isolates"

    # A second canSNP of B that lacks a base in one genome is no clean site,
    # but it is part of the set, and the finding says so.
    scheme.snp_set = [SchemeSnp("B", "", "", "G", 4, "G", "A", rank=1),
                      SchemeSnp("B", "", "", "G", 8, "G", "A", rank=2)]
    for name in clusters["B"]:
        alignment[name] = "AAAGAAAG"
    alignment[clusters["B"][0]] = "AAAGAAAN"
    findings = check_against_alignment(scheme, alignment, reference, clusters)
    single = {f.node: f.detail for f in findings if f.category == "single-site"}
    assert "The node's 1 other canSNP lack" in single["B"]
    assert "nothing else" not in single["B"]


def test_reference_derived_names_each_node_once():
    scheme = scheme_of([["R"], ["R", "A"]], [("A", "G", 2, "G", "A")])
    scheme.snp_set = [SchemeSnp("A", "", "", "G", 2, "G", "A", rank=1),
                      SchemeSnp("A", "", "", "G", 4, "G", "A", rank=2)]
    reference = Reference.from_records([("c1", "AGAGAA")])
    (found,) = [f for f in check_against_reference(scheme, reference)
                if f.category == "reference-derived"]
    assert "2 SNPs at 1 nodes" in found.detail and "branches: A." in found.detail


# -- the Markdown report -----------------------------------------------------


def test_the_markdown_report_holds_the_run_the_findings_and_the_typing():
    extra = [SchemeSnp("A", "", "", "G", 7, "G", "A", rank=2)]
    alignment = {"ref": "AAAAAAAA", "g": "AGANAGGA", "h": "AGAGAGAA"}
    scheme = chain_scheme(extra)
    _genomes, report = typed(scheme, alignment)
    findings = [problem(MODIFY, "homoplasy", "A|x", "carried outside"),
                *[recommendation(MODIFY, "spacing", f"n{i}", "close") for i in range(12)]]
    text = markdown(scheme, findings, report, [("alignment", "a.fasta")])
    assert text.startswith("# Evaluation of the 't' scheme\n")
    for heading in ("## The run", "## The scheme", "## Findings", "### Problems",
                    "## Typing of the genomes", "### Where the single rule stopped",
                    "### The canSNPs by rank", "### Nodes these genomes do not test",
                    "### Genomes the two rules type differently, or not at all"):
        assert heading in text, heading
    assert "| problem | homoplasy | modify | 1 |" in text
    assert "`A|x`: carried outside" in text
    assert "- and 2 more, in the findings table" in text
    assert "| A | missing | 1 |" in text
    for line in text.splitlines():
        if line.startswith("|"):
            assert line.endswith("|")
    assert text.isascii()


def test_without_an_alignment_the_report_has_no_typing_section():
    scheme = chain_scheme()
    text = markdown(scheme, [])
    assert "Nothing was found" in text and "## Typing" not in text


def test_audit_writes_the_report(written, dataset, tmp_path):
    outdir, _scores = written
    path = changed(dataset, tmp_path, {"s1": {PLANTED["A"][0]: "N"}})
    out = tmp_path / "report.md"
    status = cli_run("audit", "--scheme-dir", outdir, "--reference", dataset.reference,
                     "--alignment", path, "--genome-map", dataset.genome_map,
                     "--report", str(out))
    assert status in (0, 2)
    text = out.read_text()
    assert "## Typing of the genomes" in text and "## Findings" in text
    assert f"| alignment | {path} |" in text
    assert "| A | missing | 1 |" in text, "s1 stopped at A under the single rule"
    assert not [name for name in os.listdir(tmp_path) if name.startswith("typing")], (
        "the tables are written with --typing-report only")


def test_the_report_is_written_without_an_alignment_too(written, dataset, tmp_path):
    outdir, _scores = written
    out = tmp_path / "report.md"
    assert cli_run("audit", "--scheme-dir", outdir, "--reference", dataset.reference,
                   "--report", str(out)) in (0, 2)
    assert "## Typing" not in out.read_text()


def set_scheme(positions):
    """A node A whose ranked canSNPs lie at the given positions of genome G."""
    scheme = scheme_of([["R"], ["R", "A"]], [("A", "G", positions[0], "T", "G")])
    scheme.snp_set = [SchemeSnp("A", "", "", "G", position, "T", "G", rank=rank)
                      for rank, position in enumerate(positions, start=1)]
    return scheme


def test_a_node_whose_canSNPs_lie_in_one_region_is_reported():
    (found,) = check_set_region(set_scheme([127613, 3513465, 129091, 131813, 133619]), 10000)
    assert found.category == "clustered-set" and found.action == MODIFY and found.node == "A"
    assert "4 of the node's 5 canSNPs lie within 6006 bases" in found.detail
    assert "127613-133619" in found.detail


def test_canSNPs_spread_over_the_genome_or_too_few_are_not_reported():
    assert check_set_region(set_scheme([1000, 60000, 120000, 180000, 240000]), 10000) == []
    # Half of them close together is not more than half.
    assert check_set_region(set_scheme([1000, 2000, 50000, 90000]), 10000) == []
    # Two canSNPs are too few to call a set clustered.
    assert check_set_region(set_scheme([1000, 1500]), 10000) == []
    # 0 turns the check off.
    assert check_set_region(set_scheme([1000, 1100, 1200]), 0) == []
