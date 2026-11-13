"""Defects found in the fourth audit round, each pinned by the input that showed it."""

import logging

from conftest import cli_run, vcf_records
from test_review_fixes import scheme_of
from test_typing_report import run_report, snp, table

from cansnpscorer.cansnp_audit import check_callability
from cansnpscorer.cansnp_evaluator import locus_at
from cansnpscorer.helpers.data_utils import ranking
from cansnpscorer.helpers.file_utils import read_alignment, read_gff3
from cansnpscorer.hierarchy import Hierarchy
from cansnpscorer.model import Reference
from cansnpscorer.scheme import SchemeSnp
from cansnpscorer.tree_tools import place
from cansnpscorer.typing_report import (
    ANCESTRAL,
    DERIVED,
    MISSING,
    derived_off,
    read_observed,
    set_path,
    type_alignment,
)

LINEAGE = {"r": ["r"], "A": ["A", "r"], "A1": ["A1", "A", "r"], "A11": ["A11", "A1", "A", "r"],
           "B": ["B", "r"], "B1": ["B1", "B", "r"], "B11": ["B11", "B1", "B", "r"],
           "B111": ["B111", "B11", "B1", "B", "r"]}


# -- the set rule ------------------------------------------------------------


def test_one_derived_node_below_unknown_ones_does_not_outweigh_a_supported_path():
    states = {"r": DERIVED, "A": DERIVED, "A1": DERIVED, "A11": DERIVED,
              "B": MISSING, "B1": MISSING, "B11": MISSING, "B111": DERIVED}
    call, path = set_path(states, LINEAGE)
    assert call == "A11" and "B111" not in path
    assert derived_off(states, LINEAGE, call, path) == 1, "the stray node, not the path"


def test_an_unknown_node_is_still_passed_when_nothing_competes():
    call, path = set_path({"r": DERIVED, "A": MISSING, "A1": DERIVED}, LINEAGE)
    assert call == "A1" and path == {"r", "A", "A1"}


def test_two_branches_of_equal_support_give_no_call_but_do_count_as_conflict():
    states = {"r": DERIVED, "A": DERIVED, "A1": DERIVED, "B": DERIVED, "B1": DERIVED}
    call, path = set_path(states, LINEAGE)
    assert call is None and path == {"r"}
    # One of the two branches may be the genome's own: two nodes contradict it.
    assert derived_off(states, LINEAGE, call, path) == 2


def test_an_ancestral_node_still_closes_the_branch_below_it():
    states = {"r": DERIVED, "A": ANCESTRAL, "A1": DERIVED, "A11": DERIVED, "B": DERIVED}
    assert set_path(states, LINEAGE)[0] == "B"


# -- the single rule and its labels -----------------------------------------


def typed(scheme, alignment, **arguments):
    reference = Reference.from_records([("c1", alignment["ref"])])
    report = type_alignment(scheme, alignment, reference, **arguments)
    return {result.genome: result for result in report.genomes}, report


def chain_scheme(extra=()):
    """root > A > A.1, a SNP each at positions 2, 4 and 6 (G derived, A ancestral)."""
    scheme = scheme_of([["root"], ["root", "A"], ["root", "A", "A.1"], ["root", "B"]],
                       [("root", "G", 2, "G", "A"), ("A", "G", 4, "G", "A"),
                        ("A.1", "G", 6, "G", "A"), ("B", "G", 8, "G", "A")])
    scheme.snp_set = list(extra)
    return scheme


def test_a_single_call_on_another_branch_is_labelled_as_such():
    # First canSNPs: derived down root, A, A.1. B's set of three is derived by
    # majority and A's set is ancestral by majority, so the set rule calls B.
    extra = [SchemeSnp("A", "", "", "G", 10, "G", "A", rank=2),
             SchemeSnp("A", "", "", "G", 11, "G", "A", rank=3),
             SchemeSnp("B", "", "", "G", 12, "G", "A", rank=2),
             SchemeSnp("B", "", "", "G", 13, "G", "A", rank=3)]
    alignment = {"ref": "AAAAAAAAAAAAAA",
                 #        .G.G.G.A.AAGG
                 "g": "AGAGAGAAAAAGGA"}
    genomes, report = typed(chain_scheme(extra), alignment)
    result = genomes["g"]
    assert (result.single_call, result.set_call) == ("A.1", "B")
    assert result.single_stop == ("A.1", "other-branch")
    stops = {row[0]: row[-1] for row in report.node_rows}
    assert set(stops.values()) == {"0"}, "no node lost a genome of the single rule"
    summary = {metric: single for metric, single, _set in report.summary}
    assert summary["single call on another branch than the set call"] == "1"
    assert summary["single call above the set call"] == "0"


def test_the_single_rule_does_not_read_rank_2_where_rank_1_is_absent():
    scheme = chain_scheme([SchemeSnp("A", "", "", "G", 4, "G", "A", rank=2)])
    scheme.snps = [row for row in scheme.snps if row.node != "A"]
    genomes, _report = typed(scheme, {"ref": "AAAAAAAA", "g": "AGAGAGAA"})
    assert genomes["g"].single_call is None, "CanSNPer2 has no SNP at A"
    assert genomes["g"].set_call == "A.1"
    assert genomes["g"].single_stop == ("A", "missing")


def test_the_mean_depth_is_over_the_calls_made():
    _genomes, report = typed(chain_scheme(), {"ref": "AAAAAAAA", "g": "AGAGAGAA",
                                              "none": "AAAAAAAA"})
    summary = {metric: (single, by_set) for metric, single, by_set in report.summary}
    assert summary["mean depth of the calls made"] == ("3.00", "3.00")


def test_the_findings_are_written_when_the_genomes_cannot_be_typed(written, dataset,
                                                                    tmp_path, caplog):
    # The first record of the alignment is not the reference: the audit says
    # so in a finding, which a failed typing report must not take away.
    outdir, _scores = written
    records = list(read_alignment(dataset.alignment).items())
    from conftest import write_fasta

    path = str(tmp_path / "reordered.fasta")
    write_fasta(path, records[1:] + records[:1])
    output = tmp_path / "findings.tsv"
    with caplog.at_level(logging.ERROR):
        status = cli_run("audit", "--scheme-dir", outdir, "--reference", dataset.reference,
                         "--alignment", path, "--typing-report", str(tmp_path / "t"),
                         "-o", str(output))
    assert status in (1, 2) and output.exists()
    assert "No typing report was written" in caplog.text


def test_nodes_below_a_node_without_a_snp_are_not_callable():
    scheme = scheme_of(
        [["root"], ["root", "A"], ["root", "A", "A.1"], ["root", "A", "A.1", "A.1.1"]],
        [("root", "G", 2, "G", "A"), ("A.1", "G", 6, "G", "A"), ("A.1.1", "G", 8, "G", "A")])
    findings = check_callability(scheme, 3)
    unreachable = [finding for finding in findings if finding.category == "unreachable"]
    assert [finding.node for finding in unreachable] == ["A"]
    assert "A.1, A.1.1" in unreachable[0].detail
    assert any(finding.category == "min-required-hits" and finding.severity == "problem"
               for finding in findings), "no leaf can be called"


# -- evaluate ----------------------------------------------------------------


def observed_table(tmp_path, rows):
    header = ("node\trank\tgenome\tposition\tderived_base\tancestral_base\tgenomes\tcalled\t"
              "missing\tderived\tancestral\tother\tdiscordant_in_set\tancestral_on_path\t"
              "derived_off_path")
    path = tmp_path / "observed_snps.tsv"
    path.write_text(header + "\n" + "".join("\t".join(map(str, row)) + "\n" for row in rows))
    return str(path)


def test_a_snp_called_in_no_genome_was_not_observed(tmp_path):
    path = observed_table(tmp_path, [("A", 1, "T", 461, "C", "G", 6, 0, 6, 0, 0, 0, 0, 0, 0)])
    seen = read_observed(path)[("A", 461)]
    assert seen.genomes == 0 and seen.missing == 6
    never_called = snp(461, total=180, observed_genomes=seen.genomes,
                       observed_missing=seen.missing)
    untried = snp(201, total=350)
    assert ranking(untried) < ranking(never_called)


def test_two_rows_for_one_snp_keep_the_worse(tmp_path):
    path = observed_table(tmp_path, [
        ("A", 1, "T", 461, "C", "G", 6, 6, 0, 3, 3, 0, 5, 0, 0),
        ("A", 2, "T", 461, "C", "G", 6, 6, 0, 3, 3, 0, 0, 0, 0)])
    assert read_observed(path)[("A", 461)].discordant == 5


def test_a_row_with_other_alleles_is_not_taken_for_the_snp(written, dataset, tmp_path,
                                                            caplog):
    _status, prefix = run_report(written, dataset, tmp_path)
    rows = table(prefix + "_snps.tsv")
    swapped = observed_table(tmp_path, [
        [row["node"], row["rank"], "OTHER", row["position"], "N", "N", 6, 6, 0, 3, 3, 0, 4, 0,
         0] for row in rows])
    scores = str(tmp_path / "observed.vcf")
    with caplog.at_level(logging.WARNING):
        assert cli_run("evaluate", "--snps", str(tmp_path / "snps.vcf"), "--reference",
                       dataset.reference, "--genes", dataset.gff, "--repeats",
                       dataset.repeats, "-o", scores, "--observed", swapped) == 0
    assert "other alleles" in caplog.text
    assert all(record.info["OBS_DISCORDANT"] == "0" for record in vcf_records(scores))


def test_abutting_genes_are_two_loci(tmp_path):
    gff = tmp_path / "genes.gff3"
    gff.write_text("##gff-version 3\n"
                   "c1\tt\tgene\t151\t450\t.\t+\t.\tID=g1\n"
                   "c1\tt\tCDS\t151\t450\t.\t+\t0\tID=c1\n"
                   "c1\tt\tgene\t451\t650\t.\t-\t.\tID=g2\n"
                   "c1\tt\tgene\t640\t700\t.\t+\t.\tID=g3\n")
    model = read_gff3(str(gff))
    assert model.genes["c1"] == [(150, 699)], "joined, for the gene term"
    loci = model.loci["c1"]
    assert loci == [(150, 449), (450, 649), (639, 699)]
    starts, longest = [start for start, _end in loci], 300

    def at(position):
        return locus_at(position - 1, loci, starts, longest)

    assert at(450) == (150, 449) and at(451) == (450, 649)
    assert at(645) == (450, 649), "of two overlapping genes, the first"
    assert at(700) == (639, 699) and at(701) is None and at(100) is None


def test_a_cutoff_no_score_can_pass_is_said(written, dataset, tmp_path, caplog):
    with caplog.at_level(logging.WARNING):
        assert cli_run("evaluate", "--snps", str(tmp_path / "snps.vcf"), "--reference",
                       dataset.reference, "--genes", dataset.gff, "--repeats",
                       dataset.repeats, "-o", str(tmp_path / "s.vcf"), "--profile", "assay",
                       "--score-weight", "gc_score=0") == 0
    assert "No canSNP can pass" in caplog.text


def test_a_copy_number_that_could_not_be_counted_ranks_after_unique_sequence():
    assert min([snp(1, copies=2), snp(2, copies=0), snp(3, copies=1)],
               key=ranking).copies == 1
    assert ranking(snp(1, total=100, copies=1)) < ranking(snp(2, total=600, copies=0))


# -- place -------------------------------------------------------------------


def test_genomes_without_a_base_at_the_parents_snp_do_not_break_the_order():
    scheme = scheme_of([["root"], ["root", "X"], ["root", "X", "X.1"]],
                       [("X", "G", 3, "T", "G"), ("X.1", "G", 6, "A", "C")])
    reference = Reference.from_records([("c1", "ACGTACGT")])
    alignment = {"r": "ACGTACGT", "x1": "ACTTACGT", "x2": "ACTTACGT", "x3": "ACTTACGT",
                 "x4": "AC-TAAGT", "x5": "AC-TAAGT", "y1": "ACGTACGT", "y2": "ACGTACGT"}
    tree = Hierarchy.from_newick("(r,(((x1,x2,x3),(x4,x5)),(y1,y2)));", support_labels=True)
    placement = place(scheme, tree, alignment, reference)
    assert not [finding for finding in placement.findings
                if finding.category == "placement-order"]
    assert sorted(placement.clusters["X.1"]) == ["x4", "x5"]


def test_one_outgroup_of_an_unrooted_tree_stays_outside_the_root():
    scheme = scheme_of([["root"]], [("root", "G", 3, "T", "G")])
    reference = Reference.from_records([("c1", "ACGTACGT")])
    alignment = {"r": "ACTTACGT", "a": "ACTTACGT", "b": "ACTTACGT", "c": "ACTTACGT",
                 "o": "ACGTACGT"}
    alignment = {"ref": "ACGTACGT", **alignment}
    tree = Hierarchy.from_newick("(r,a,(o,(b,c)));", support_labels=True)
    placement = place(scheme, tree, alignment, reference)
    assert sorted(placement.clusters["root"]) == ["a", "b", "c", "r"]


def test_a_node_is_placed_by_its_best_cansnp_not_only_its_first():
    scheme = scheme_of([["root"], ["root", "X"]], [("X", "G", 3, "T", "G")])
    scheme.snp_set = [SchemeSnp("X", "", "", "G", 6, "A", "C", rank=2)]
    reference = Reference.from_records([("c1", "ACGTACGT")])
    # Rank 1 is derived in x1 and in y2; rank 2 in exactly x1, x2 and x3.
    alignment = {"r": "ACGTACGT", "x1": "ACTTAAGT", "x2": "ACGTAAGT", "x3": "ACGTAAGT",
                 "y1": "ACGTACGT", "y2": "ACTTACGT"}
    tree = Hierarchy.from_newick("(r,((x1,x2,x3),(y1,y2)));", support_labels=True)
    placement = place(scheme, tree, alignment, reference)
    assert sorted(placement.clusters["X"]) == ["x1", "x2", "x3"]
