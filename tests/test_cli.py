"""The command line, run as a user runs it.

These go through cansnp_scorer.main, so they cover argument handling, the exit
status of each subcommand and the files they leave behind.
"""

import csv
import os

import pytest
from conftest import PLANTED

from cansnpscorer.cansnp_scorer import main, score_weights
from cansnpscorer.helpers.scoring_methods import DEFAULT_WEIGHTS


def run(monkeypatch, *arguments):
    monkeypatch.setattr("sys.argv", ["CanSNPCreator", *arguments])
    return main()


def extract(monkeypatch, dataset, output, hierarchy=None, groups=None):
    return run(
        monkeypatch, "extract",
        "--genome-map", groups or dataset.genome_map,
        "--hierarchy", hierarchy or dataset.hierarchy,
        "--alignment", dataset.alignment,
        "--reference", dataset.reference,
        "-o", output,
    )


def evaluate(monkeypatch, dataset, snp_file, output, *extra):
    return run(
        monkeypatch, "evaluate",
        "--snps", snp_file,
        "--reference", dataset.reference,
        "--genes", dataset.gff,
        "--repeats", dataset.trf,
        "--hierarchy", dataset.hierarchy,
        "-o", output,
        "--cutoff", "0",
        *extra,
    )


# -- the whole run -----------------------------------------------------------


@pytest.mark.parametrize("which", ["table", "newick"])
def test_extract_then_evaluate(monkeypatch, dataset, tmp_path, which):
    snp_file = str(tmp_path / "snps.tsv")
    scores = str(tmp_path / "scores.csv")

    hierarchy = dataset.hierarchy if which == "table" else dataset.hierarchy_newick
    assert extract(monkeypatch, dataset, snp_file, hierarchy, dataset.genome_map) == 0
    assert evaluate(monkeypatch, dataset, snp_file, scores) == 0

    rows = list(csv.DictReader(open(scores)))
    assert {row["Group"] for row in rows} == set(PLANTED)

    # Every planted canSNP is reported on its own contig, at its 1-based
    # position, with the polarity it was planted with.
    for row in rows:
        column, ancestral, derived = PLANTED[row["Group"]]
        contig, position = dataset.contig_of(column)
        assert (row["Contig"], int(row["Position"])) == (contig, position)
        assert (row["Derived"], row["Ancestral"]) == (derived, ancestral)


def test_the_snp_file_has_a_header_and_one_row_per_snp(monkeypatch, dataset, tmp_path):
    snp_file = str(tmp_path / "snps.tsv")
    assert extract(monkeypatch, dataset, snp_file) == 0
    lines = open(snp_file).read().splitlines()
    assert lines[0] == "contig\tposition\tnode\tderived_base\tancestral_base"
    assert len(lines) == len(PLANTED) + 2, "the planted canSNPs, one incongruent SNP, a header"


def test_the_evaluation_can_be_narrowed_to_the_best_snp_per_node(monkeypatch, dataset, tmp_path):
    snp_file = str(tmp_path / "snps.tsv")
    scores = str(tmp_path / "scores.csv")
    assert extract(monkeypatch, dataset, snp_file) == 0
    assert evaluate(monkeypatch, dataset, snp_file, scores, "--top-n", "1") == 0
    assert len(list(csv.DictReader(open(scores)))) == len(PLANTED)


def test_an_existing_score_file_is_backed_up(monkeypatch, dataset, tmp_path):
    snp_file = str(tmp_path / "snps.tsv")
    scores = str(tmp_path / "scores.csv")
    assert extract(monkeypatch, dataset, snp_file) == 0
    assert evaluate(monkeypatch, dataset, snp_file, scores) == 0
    assert evaluate(monkeypatch, dataset, snp_file, scores) == 0
    backups = [name for name in os.listdir(tmp_path) if name.endswith(".backup")]
    assert len(backups) == 1


# -- failures ----------------------------------------------------------------


def test_a_missing_input_file_exits_non_zero(monkeypatch, dataset, tmp_path):
    status = run(
        monkeypatch, "extract",
        "--genome-map", dataset.genome_map,
        "--hierarchy", dataset.hierarchy,
        "--alignment", str(tmp_path / "absent.fasta"),
        "--reference", dataset.reference,
        "-o", str(tmp_path / "snps.tsv"),
    )
    assert status == 1


def test_a_reference_that_is_not_the_aligned_one_exits_two(monkeypatch, dataset, tmp_path):
    other = tmp_path / "other.fasta"
    other.write_text(">c1\n" + "A" * 1000 + "\n")
    monkeypatch.setattr("sys.argv", [
        "CanSNPCreator", "extract",
        "--genome-map", dataset.genome_map,
        "--hierarchy", dataset.hierarchy,
        "--alignment", dataset.alignment,
        "--reference", str(other),
        "-o", str(tmp_path / "snps.tsv"),
    ])
    assert main() == 2, "the inputs are inconsistent, so nothing is written"
    assert not os.path.exists(tmp_path / "snps.tsv")


def test_a_snp_file_from_another_reference_exits_non_zero(monkeypatch, dataset, tmp_path):
    # Neither allele is the base the reference carries there, which is what a
    # SNP file produced against something else looks like.
    rows = []
    for position in (1, 2):
        others = [base for base in "ACGT" if base != dataset.sequences["ref"][position - 1]]
        rows.append(f"c1\t{position}\tA.1\t{others[0]}\t{others[1]}")
    snp_file = tmp_path / "snps.tsv"
    snp_file.write_text(
        "contig\tposition\tnode\tderived_base\tancestral_base\n" + "\n".join(rows) + "\n")
    assert evaluate(monkeypatch, dataset, str(snp_file), str(tmp_path / "scores.csv")) == 1


def test_no_subcommand_prints_help_and_exits_two(monkeypatch, capsys):
    assert run(monkeypatch) == 2
    assert "extract" in capsys.readouterr().out


# -- arguments ---------------------------------------------------------------


def test_score_weights_default_to_the_doubled_gene_term():
    assert score_weights(None) == DEFAULT_WEIGHTS


def test_score_weights_are_parsed():
    assert score_weights(["gene_score=1", "gc_score=0.5"]) == {"gene_score": 1.0, "gc_score": 0.5}


def test_an_unknown_score_term_is_refused():
    with pytest.raises(Exception, match="unknown score term"):
        score_weights(["genes=2"])


def test_a_non_numeric_weight_is_refused():
    with pytest.raises(Exception, match="must be a number"):
        score_weights(["gene_score=high"])


def test_weighting_the_gene_term_once_lowers_the_total(monkeypatch, dataset, tmp_path):
    snp_file = str(tmp_path / "snps.tsv")
    doubled = str(tmp_path / "doubled.csv")
    plain = str(tmp_path / "plain.csv")
    assert extract(monkeypatch, dataset, snp_file) == 0
    assert evaluate(monkeypatch, dataset, snp_file, doubled) == 0
    assert evaluate(monkeypatch, dataset, snp_file, plain, "--score-weight", "gene_score=1") == 0

    def totals(path):
        return {row["Group"]: int(row["Score"]) for row in csv.DictReader(open(path))}

    for node, total in totals(plain).items():
        assert total <= totals(doubled)[node]
    assert totals(plain) != totals(doubled)
