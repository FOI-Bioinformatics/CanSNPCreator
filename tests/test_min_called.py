"""extract --min-called: SNPs in columns where some genomes have no base."""

from types import SimpleNamespace

from conftest import GENOMES, PLANTED, vcf_records, write_fasta
from test_cli import evaluate, run

from cansnpscorer.helpers.data_utils import (
    called_needed,
    canSNP_extractor,
    ranking,
    select_per_node,
    snp_extractor,
)

ALIGNMENT = {
    "ref": "AAAC",
    "g1": "TTAC",
    "g2": "TNAC",
    "g3": "AAA-",
    "g4": "AAAC",
}
GROUPS = {"X": ["g1", "g2"], "Y": ["g3"]}


def test_called_needed_rounds_up():
    assert called_needed(10, 0.9) == 9
    assert called_needed(5, 0.9) == 5
    assert called_needed(2, 0.5) == 1
    assert called_needed(10, 1.0) == 10


def test_default_leaves_out_any_column_with_a_missing_base():
    snps = snp_extractor(ALIGNMENT)
    assert [snp["index"] for snp in snps] == [0]
    assert snps[0]["missing"] == 0


def test_a_fraction_keeps_a_column_with_enough_genomes_called():
    snps = snp_extractor(ALIGNMENT, min_called=0.8)
    assert [snp["index"] for snp in snps] == [0, 1]
    assert snps[1]["variants"] == {"T": 1, "A": 3}
    assert snps[1]["missing"] == 1
    # Too few called for a stricter fraction.
    assert [snp["index"] for snp in snp_extractor(ALIGNMENT, min_called=0.9)] == [0]


def test_a_column_with_one_base_among_the_called_genomes_is_no_snp():
    # Column 3 holds C and a gap, so only one base among the called genomes.
    assert 3 not in [snp["index"] for snp in snp_extractor(ALIGNMENT, min_called=0.5)]


def test_the_called_genomes_of_a_node_decide_it():
    snps = snp_extractor(ALIGNMENT, min_called=0.5)
    found = {snp.column: snp for snp in canSNP_extractor(ALIGNMENT, GROUPS, snps, 0.5)}
    assert found[1].node == "X" and found[1].derived == "T" and found[1].ancestral == "A"
    assert found[1].missing == 1


def test_a_node_with_too_few_called_genomes_is_not_marked():
    snps = snp_extractor(ALIGNMENT, min_called=0.8)
    # X has two genomes and one is missing; 0.8 of two is two.
    found = {snp.column: snp for snp in canSNP_extractor(ALIGNMENT, GROUPS, snps, 0.8)}
    assert not found[1].is_canonical


def test_a_derived_base_outside_the_node_still_breaks_it():
    alignment = dict(ALIGNMENT, g4="ATAC")
    snps = snp_extractor(alignment, min_called=0.5)
    found = {snp.column: snp for snp in canSNP_extractor(alignment, GROUPS, snps, 0.5)}
    assert not found[1].is_canonical


def test_ranking_takes_the_complete_snp_before_a_better_scoring_incomplete_one():
    complete = SimpleNamespace(node="X", total=50, contig="c", position=10, missing=0)
    partial = SimpleNamespace(node="X", total=90, contig="c", position=5000, missing=2)
    assert ranking(complete) < ranking(partial)
    assert select_per_node([partial, complete], top=1)["X"] == [complete]


def _with_missing_base(dataset, tmp_path):
    """The dataset's alignment with s2's base at the root's canSNP replaced by N."""
    column = PLANTED["root"][0]
    sequences = dict(dataset.sequences)
    sequences["s2"] = sequences["s2"][:column] + "N" + sequences["s2"][column + 1:]
    path = str(tmp_path / "missing.fasta")
    write_fasta(path, [(name, sequences[name]) for name in GENOMES])
    return path


def _extract(monkeypatch, dataset, alignment, output, *extra):
    return run(monkeypatch, "extract", "--genome-map", dataset.genome_map,
               "--hierarchy", dataset.hierarchy, "--alignment", alignment,
               "--reference", dataset.reference, "-o", output, *extra)


def test_cli_min_called_recovers_a_canSNP_with_a_missing_base(monkeypatch, dataset, tmp_path):
    alignment = _with_missing_base(dataset, tmp_path)
    strict = str(tmp_path / "strict.vcf")
    assert _extract(monkeypatch, dataset, alignment, strict) == 0
    assert "root" not in {record.node for record in vcf_records(strict)}

    relaxed = str(tmp_path / "relaxed.vcf")
    assert _extract(monkeypatch, dataset, alignment, relaxed, "--min-called", "0.8") == 0
    nodes = {record.node: record for record in vcf_records(relaxed) if record.node}
    assert set(nodes) == set(PLANTED)
    assert nodes["root"].info["MISSING"] == "1"
    assert nodes["A"].info["MISSING"] == "0"

    scores = str(tmp_path / "scores.vcf")
    assert evaluate(monkeypatch, dataset, relaxed, scores) == 0
    scored = {record.node: record for record in vcf_records(scores)}
    assert scored["root"].info["MISSING"] == "1"


def test_cli_refuses_a_fraction_out_of_range(monkeypatch, dataset, tmp_path, capsys):
    for value in ("0", "1.5"):
        try:
            _extract(monkeypatch, dataset, dataset.alignment, str(tmp_path / "x.vcf"),
                     "--min-called", value)
        except SystemExit as stop:
            assert stop.code == 2
        else:
            raise AssertionError(f"--min-called {value} was accepted")


def test_a_node_is_not_marked_when_the_isolates_above_it_all_lack_a_base():
    # P holds X (x1-x3) and Y (y1, y2). A SNP of P with both Y genomes missing
    # looks like a SNP of X, and P itself has too few genomes called to match,
    # so nothing else flags the column.
    alignment = {"out": "AA", "x1": "TT", "x2": "TT", "x3": "TT", "y1": "NA", "y2": "NA",
                 **{f"o{i}": "AA" for i in range(14)}}
    groups = {"P": ["x1", "x2", "x3", "y1", "y2"], "X": ["x1", "x2", "x3"], "Y": ["y1", "y2"]}
    snps = snp_extractor(alignment, min_called=0.75)
    found = {snp.column: snp for snp in canSNP_extractor(alignment, groups, snps, 0.75)}
    assert found[0].node is None, "X or P: the column cannot tell"
    assert found[1].node == "X", "with Y called ancestral it is X's"
