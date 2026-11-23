"""The command line, run as a user runs it.

These go through cansnp_scorer.main, so they cover argument handling, the exit
status of each subcommand and the files they leave behind.
"""

import os
from pathlib import Path

import pytest
from conftest import PLANTED, vcf_lines, vcf_records

from cansnpscorer.cansnp_scorer import main, score_weights
from cansnpscorer.helpers.scoring_methods import DEFAULT_WEIGHTS, TYPING_WEIGHTS


def run(monkeypatch, *arguments):
    monkeypatch.setattr("sys.argv", ["CanSNPCreator", *arguments])
    return main()


def extract(monkeypatch, dataset, output, hierarchy=None, groups=None, *extra):
    return run(
        monkeypatch,
        "extract",
        "--genome-map",
        groups or dataset.genome_map,
        "--hierarchy",
        hierarchy or dataset.hierarchy,
        "--alignment",
        dataset.alignment,
        "--reference",
        dataset.reference,
        "-o",
        output,
        *extra,
    )


def placed_everywhere(dataset, tmp_path):
    """The planted genome map with the reference placed at A.1 too, so that no
    genome of the alignment sits outside every node."""
    path = tmp_path / "all_placed.tsv"
    lines = Path(dataset.genome_map).read_text().splitlines()
    path.write_text("\n".join([*lines, "ref\tA.1"]) + "\n")
    return str(path)


def evaluate(monkeypatch, dataset, snp_file, output, *extra):
    return run(
        monkeypatch,
        "evaluate",
        "--snps",
        snp_file,
        "--reference",
        dataset.reference,
        "--genes",
        dataset.gff,
        "--repeats",
        dataset.repeats,
        "--hierarchy",
        dataset.hierarchy,
        "-o",
        output,
        "--cutoff",
        "0",
        *extra,
    )


# -- the whole run -----------------------------------------------------------


@pytest.mark.parametrize("which", ["table", "newick"])
def test_extract_then_evaluate(monkeypatch, dataset, tmp_path, which):
    snp_file = str(tmp_path / "snps.vcf")
    scores = str(tmp_path / "scores.vcf")

    hierarchy = dataset.hierarchy if which == "table" else dataset.hierarchy_newick
    assert extract(monkeypatch, dataset, snp_file, hierarchy, dataset.genome_map) == 0
    assert evaluate(monkeypatch, dataset, snp_file, scores) == 0

    records = vcf_records(scores)
    assert {record.node for record in records} == set(PLANTED)

    # Every planted canSNP is reported on its own contig, at its 1-based
    # position, with the polarity it was planted with, and REF is the
    # reference's base.
    for record in records:
        column, ancestral, derived = PLANTED[record.node]
        contig, position = dataset.contig_of(column)
        assert (record.chrom, record.pos) == (contig, position)
        assert (record.derived, record.ancestral) == (derived, ancestral)
        assert record.ref == dataset.sequences["ref"][column]
        assert record.passed


def test_each_snp_records_the_nearest_variation_among_the_placed_genomes(monkeypatch, dataset, tmp_path):
    snp_file = str(tmp_path / "snps.vcf")
    assert extract(monkeypatch, dataset, snp_file) == 0
    header, _records = vcf_lines(snp_file)
    assert any(line.startswith("##INFO=<ID=NEAREST,") for line in header)
    # The planted sites lie 100 columns apart, the backbone is the same in
    # every genome, and the root's site does not vary among the placed
    # genomes, so each canSNP's nearest variation is 100 bases away; B.2 sits
    # at the start of c2, whose nearest is the incongruent site at 800.
    nearest = {record.node: record.info["NEAREST"] for record in vcf_records(snp_file) if record.node}
    assert nearest == dict.fromkeys(PLANTED, "100")


def test_the_snps_are_written_as_vcf(monkeypatch, dataset, tmp_path):
    snp_file = str(tmp_path / "snps.vcf")
    assert extract(monkeypatch, dataset, snp_file) == 0
    header, records = vcf_lines(snp_file)
    assert header[0] == "##fileformat=VCFv4.3"
    assert "##contig=<ID=c1,length=700>" in header and "##contig=<ID=c2,length=300>" in header
    assert any(line.startswith("##INFO=<ID=AA,") for line in header)
    assert header[-1] == "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    assert len(records) == len(PLANTED) + 1, "the planted canSNPs and one incongruent SNP"
    nodes = [record.node for record in vcf_records(snp_file)]
    assert nodes.count(None) == 1, "the incongruent SNP marks no node and carries no AA"


def test_the_evaluation_can_be_narrowed_to_the_best_snp_per_node(monkeypatch, dataset, tmp_path):
    snp_file = str(tmp_path / "snps.vcf")
    scores = str(tmp_path / "scores.vcf")
    assert extract(monkeypatch, dataset, snp_file) == 0
    assert evaluate(monkeypatch, dataset, snp_file, scores, "--top-n", "1") == 0
    assert sum(record.passed for record in vcf_records(scores)) == len(PLANTED)


def test_failing_canSNPs_are_kept_with_a_filter(monkeypatch, dataset, tmp_path):
    snp_file = str(tmp_path / "snps.vcf")
    scores = str(tmp_path / "scores.vcf")
    assert extract(monkeypatch, dataset, snp_file) == 0
    assert evaluate(monkeypatch, dataset, snp_file, scores, "--cutoff", "600") == 0
    records = vcf_records(scores)
    assert len(records) == len(PLANTED), "nothing is dropped"
    assert all(record.filters == ["LowScore"] for record in records)
    header, _records = vcf_lines(scores)
    assert any("cutoff of 600" in line for line in header if line.startswith("##FILTER"))


def test_an_existing_output_is_replaced_without_a_backup(monkeypatch, dataset, tmp_path):
    # Every subcommand overwrites its output, as command line tools do.
    snp_file = str(tmp_path / "snps.vcf")
    scores = str(tmp_path / "scores.vcf")
    assert extract(monkeypatch, dataset, snp_file) == 0
    assert evaluate(monkeypatch, dataset, snp_file, scores) == 0
    assert evaluate(monkeypatch, dataset, snp_file, scores) == 0
    assert not [name for name in os.listdir(tmp_path) if name.endswith(".backup")]


# -- failures ----------------------------------------------------------------


def test_a_missing_input_file_exits_non_zero(monkeypatch, dataset, tmp_path):
    status = run(
        monkeypatch,
        "extract",
        "--genome-map",
        dataset.genome_map,
        "--hierarchy",
        dataset.hierarchy,
        "--alignment",
        str(tmp_path / "absent.fasta"),
        "--reference",
        dataset.reference,
        "-o",
        str(tmp_path / "snps.vcf"),
    )
    assert status == 1


def test_a_reference_that_is_not_the_aligned_one_exits_two(monkeypatch, dataset, tmp_path):
    other = tmp_path / "other.fasta"
    other.write_text(">c1\n" + "A" * 1000 + "\n")
    monkeypatch.setattr(
        "sys.argv",
        [
            "CanSNPCreator",
            "extract",
            "--genome-map",
            dataset.genome_map,
            "--hierarchy",
            dataset.hierarchy,
            "--alignment",
            dataset.alignment,
            "--reference",
            str(other),
            "-o",
            str(tmp_path / "snps.vcf"),
        ],
    )
    assert main() == 2, "the inputs are inconsistent, so nothing is written"
    assert not os.path.exists(tmp_path / "snps.vcf")


def test_snps_from_another_reference_are_refused(monkeypatch, dataset, tmp_path):
    # The ##contig lines record the reference the SNPs were found against.
    snp_file = str(tmp_path / "snps.vcf")
    assert extract(monkeypatch, dataset, snp_file) == 0
    text = Path(snp_file).read_text().replace("##contig=<ID=c2,length=300>", "##contig=<ID=c2,length=301>")
    other = tmp_path / "other.vcf"
    other.write_text(text)
    assert evaluate(monkeypatch, dataset, str(other), str(tmp_path / "scores.vcf")) == 1


def test_a_file_that_is_not_vcf_is_refused(monkeypatch, dataset, tmp_path):
    snp_file = tmp_path / "snps.tsv"
    snp_file.write_text("contig\tposition\tnode\tderived_base\tancestral_base\n")
    assert evaluate(monkeypatch, dataset, str(snp_file), str(tmp_path / "scores.vcf")) == 1


def test_no_subcommand_prints_help_and_exits_two(monkeypatch, capsys):
    assert run(monkeypatch) == 2
    assert "extract" in capsys.readouterr().out


# -- arguments ---------------------------------------------------------------


def test_score_weights_default_to_those_of_the_typing_profile():
    assert score_weights(None) == TYPING_WEIGHTS


def test_the_assay_profile_doubles_the_gene_term():
    assert score_weights(None, "assay") == DEFAULT_WEIGHTS
    assert score_weights(["gc_score=0.5"], "assay") == {"gene_score": 2.0, "gc_score": 0.5}


def test_score_weights_are_parsed():
    assert score_weights(["gene_score=1", "gc_score=0.5"]) == {"gene_score": 1.0, "gc_score": 0.5}


def test_an_unknown_score_term_is_refused():
    with pytest.raises(Exception, match="unknown score term"):
        score_weights(["genes=2"])


def test_a_non_numeric_weight_is_refused():
    with pytest.raises(Exception, match="must be a number"):
        score_weights(["gene_score=high"])


def test_weighting_the_gene_term_once_lowers_the_total(monkeypatch, dataset, tmp_path):
    snp_file = str(tmp_path / "snps.vcf")
    doubled = str(tmp_path / "doubled.vcf")
    plain = str(tmp_path / "plain.vcf")
    assert extract(monkeypatch, dataset, snp_file) == 0
    assert evaluate(monkeypatch, dataset, snp_file, doubled, "--profile", "assay") == 0
    assert (
        evaluate(
            monkeypatch, dataset, snp_file, plain, "--profile", "assay", "--score-weight", "gene_score=1"
        )
        == 0
    )

    def totals(path):
        return {record.node: int(record.info["SCORE"]) for record in vcf_records(path)}

    for node, total in totals(plain).items():
        assert total <= totals(doubled)[node]
    assert totals(plain) != totals(doubled)


def test_a_canSNP_in_a_masked_region_is_kept_with_a_filter(monkeypatch, dataset, tmp_path):
    # Mask the region holding A.1's canSNP, as Gubbins would report a
    # recombinant region.
    snp_file = str(tmp_path / "snps.vcf")
    scores = str(tmp_path / "scores.vcf")
    column = PLANTED["A.1"][0]
    masks = tmp_path / "recombination.bed"
    masks.write_text(f"c1\t{column - 10}\t{column + 10}\n")
    assert extract(monkeypatch, dataset, snp_file) == 0
    assert evaluate(monkeypatch, dataset, snp_file, scores, "--masks", str(masks)) == 0
    filters = {record.node: record.filters for record in vcf_records(scores)}
    assert filters["A.1"] == ["Masked"]
    assert all(not found for node, found in filters.items() if node != "A.1")


def test_masks_on_another_genome_are_refused(monkeypatch, dataset, tmp_path):
    snp_file = str(tmp_path / "snps.vcf")
    masks = tmp_path / "masks.bed"
    masks.write_text("NC_000001.1\t0\t10\n")
    assert extract(monkeypatch, dataset, snp_file) == 0
    assert evaluate(monkeypatch, dataset, snp_file, str(tmp_path / "s.vcf"), "--masks", str(masks)) == 1


def test_genomes_the_alignment_lacks_are_left_out_with_a_warning(monkeypatch, dataset, tmp_path, caplog):
    # A RepGenR2 genome map names every genome, its alignment often only the
    # representatives.
    genome_map = tmp_path / "genome_map.tsv"
    genome_map.write_text(Path(dataset.genome_map).read_text() + "GCF_000999999.1\tA.1\n")
    snp_file = str(tmp_path / "snps.vcf")
    assert extract(monkeypatch, dataset, snp_file, groups=str(genome_map)) == 0
    assert "1 genomes that the alignment does not hold" in caplog.text


# -- the outgroup of extract -----------------------------------------------------


def test_a_declared_outgroup_is_recorded_in_the_vcf(monkeypatch, dataset, tmp_path):
    output = str(tmp_path / "snps.vcf")
    assert extract(monkeypatch, dataset, output, None, None, "--outgroup", "ref") == 0
    header, _records = vcf_lines(output)
    assert "##outgroup=ref" in header


def test_the_inferred_outgroup_is_recorded_too(monkeypatch, dataset, tmp_path):
    output = str(tmp_path / "snps.vcf")
    assert extract(monkeypatch, dataset, output) == 0
    assert "##outgroup=ref" in vcf_lines(output)[0]


def test_no_outgroup_at_all_stops_extract(monkeypatch, dataset, tmp_path, caplog):
    output = tmp_path / "snps.vcf"
    status = extract(monkeypatch, dataset, str(output), None, placed_everywhere(dataset, tmp_path))
    assert status == 2
    assert not output.exists(), "nothing is written"
    assert "--allow-no-outgroup" in caplog.text


def test_allow_no_outgroup_takes_the_reference_base_and_says_so(monkeypatch, dataset, tmp_path):
    output = str(tmp_path / "snps.vcf")
    status = extract(
        monkeypatch, dataset, output, None, placed_everywhere(dataset, tmp_path), "--allow-no-outgroup"
    )
    assert status == 0
    assert "##outgroup=none" in vcf_lines(output)[0]


def test_an_outgroup_the_alignment_lacks_stops_extract(monkeypatch, dataset, tmp_path, caplog):
    output = tmp_path / "snps.vcf"
    status = extract(monkeypatch, dataset, str(output), None, None, "--outgroup", "GCF_000000001.1")
    assert status == 2
    assert not output.exists()
    assert "GCF_000000001.1" in caplog.text
