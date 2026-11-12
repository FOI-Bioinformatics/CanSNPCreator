"""Export of a scored scheme as CanSNPer2-data source files.

The rules checked here are the ones that repository's source/README.md sets out
and its validate.py enforces. The last test runs that validate.py itself, when
the repository is checked out beside this one.
"""

import os
import shutil
import subprocess
import sys

import pytest
from conftest import PLANTED, vcf_rewrite

from cansnpscorer.cansnp_export import best_per_node
from cansnpscorer.cansnp_export import main as export_main
from cansnpscorer.cansnp_scorer import main as cli
from cansnpscorer.helpers.scoring_methods import TERMS, ScoredSnp

METADATA = dict(
    organism="testus_organismus",
    genome="TEST1",
    strain="testus",
    citation="Sjodin2026",
    genbank_id="GCA_000007125.1",
    refseq_id="GCF_000007125.1",
    assembly_name="ASM712v1",
)


def score(node, contig="c1", position=100, total=500, derived="G", ancestral="C"):
    return ScoredSnp(contig=contig, position=position, node=node, total=total,
                     terms=dict.fromkeys(TERMS, 100), derived=derived, ancestral=ancestral)


@pytest.fixture
def scheme(dataset, tmp_path):
    """Run extract and evaluate, then return where a scheme can be written."""
    snp_file = str(tmp_path / "snps.vcf")
    scores = str(tmp_path / "scores.vcf")
    assert cli_run("extract", "--genome-map", dataset.genome_map, "--hierarchy", dataset.hierarchy,
                   "--alignment", dataset.alignment, "--reference", dataset.reference,
                   "-o", snp_file) == 0
    assert cli_run("evaluate", "--snps", snp_file, "--reference", dataset.reference,
                   "--genes", dataset.gff, "--repeats", dataset.trf,
                   "--hierarchy", dataset.hierarchy, "-o", scores,
                   "--cutoff", "0") == 0
    return scores, str(tmp_path / "scheme")


def cli_run(*arguments):
    argv = sys.argv
    try:
        sys.argv = ["CanSNPCreator", *arguments]
        return cli()
    finally:
        sys.argv = argv


def read(path):
    with open(path) as handle:
        return handle.read()


# -- choosing one SNP per node -----------------------------------------------


def test_the_highest_scoring_snp_of_each_node_is_chosen():
    scores = [score("A", position=10, total=300), score("A", position=20, total=500),
              score("B", position=30, total=400)]
    best = best_per_node(scores)
    assert (best["A"].position, best["B"].position) == (20, 30)


def test_a_tie_is_broken_by_contig_and_position():
    scores = [score("A", contig="c2", position=5, total=500),
              score("A", contig="c1", position=50, total=500)]
    assert best_per_node(scores)["A"].contig == "c1"


# -- the written files -------------------------------------------------------


def test_the_scheme_is_written_and_is_complete(dataset, scheme):
    scores, outdir = scheme
    assert export_main(scores, dataset.hierarchy, dataset.reference, outdir, **METADATA) == 0

    tree = read(os.path.join(outdir, "testus_organismus_tree.txt"))
    snp = read(os.path.join(outdir, "testus_organismus_snp.txt"))
    references = read(os.path.join(outdir, "testus_organismus_references.txt"))

    assert tree.splitlines()[0] == "root", "the first line is the root alone"
    assert all(line.startswith("root") for line in tree.splitlines())
    assert "root\tB\tB.2" in tree

    rows = [line.split("\t") for line in snp.splitlines()[1:]]
    assert {row[0] for row in rows} == set(PLANTED), "one SNP per node"
    assert len(rows) == len(tree.splitlines())

    assert references.splitlines()[1].split("\t")[0] == "TEST1"
    for text in (tree, snp, references):
        assert ";" not in text
        assert text.endswith("\n")
        assert not any(line != line.rstrip() for line in text.splitlines())


def test_positions_are_counted_over_the_concatenated_genome(dataset, scheme):
    # The SNP file has no contig column, so a SNP on the second contig is
    # offset by the length of the first, as the published schemes do it.
    scores, outdir = scheme
    assert export_main(scores, dataset.hierarchy, dataset.reference, outdir, **METADATA) == 0
    positions = {}
    for line in read(os.path.join(outdir, "testus_organismus_snp.txt")).splitlines()[1:]:
        fields = line.split("\t")
        positions[fields[0]] = int(fields[4])

    for node, (column, _ancestral, _derived) in PLANTED.items():
        assert positions[node] == column + 1, node
    assert positions["B.2"] == 701, "the first base of the second contig"


def test_the_polarity_survives_the_round_trip(dataset, scheme):
    scores, outdir = scheme
    assert export_main(scores, dataset.hierarchy, dataset.reference, outdir, **METADATA) == 0
    for line in read(os.path.join(outdir, "testus_organismus_snp.txt")).splitlines()[1:]:
        node, _strain, _reference, _genome, _position, derived, ancestral = line.split("\t")
        planted_ancestral, planted_derived = PLANTED[node][1], PLANTED[node][2]
        assert (derived, ancestral) == (planted_derived, planted_ancestral), node


# -- what is refused ---------------------------------------------------------


def test_a_node_without_a_cansnp_refuses_the_scheme(dataset, scheme, tmp_path, caplog):
    # Every node of a CanSNPer2 scheme carries exactly one SNP, so an
    # incomplete scheme is refused rather than written.
    scores, outdir = scheme
    trimmed = vcf_rewrite(scores, tmp_path / "fewer.vcf",
                          keep=lambda fields: fields[:2] != ["c1", "601"])

    assert export_main(trimmed, dataset.hierarchy, dataset.reference, outdir, **METADATA) == 2
    assert "nodes have no canSNP" in caplog.text
    assert not os.path.exists(os.path.join(outdir, "testus_organismus_tree.txt"))


def test_pruning_leaves_those_nodes_out(dataset, scheme, tmp_path):
    scores, outdir = scheme
    trimmed = vcf_rewrite(scores, tmp_path / "fewer.vcf",
                          keep=lambda fields: fields[:2] != ["c1", "201"])

    assert export_main(trimmed, dataset.hierarchy, dataset.reference, outdir,
                       prune=True, **METADATA) == 0
    tree = read(os.path.join(outdir, "testus_organismus_tree.txt"))
    assert "\tA\t" not in tree, "the node without a canSNP is gone"
    assert "root\tA.1" in tree, "its children hang from the nearest node that has one"


def test_a_root_without_a_cansnp_cannot_be_pruned(dataset, scheme, tmp_path, caplog):
    scores, outdir = scheme
    trimmed = vcf_rewrite(scores, tmp_path / "fewer.vcf",
                          keep=lambda fields: fields[:2] != ["c1", "101"])

    assert export_main(trimmed, dataset.hierarchy, dataset.reference, outdir,
                       prune=True, **METADATA) == 2
    assert "has no canSNP, and a scheme cannot be written" in caplog.text


@pytest.mark.parametrize("field,value,expected", [
    ("genome", "TEST_1", "underscore"),
    ("genbank_id", "GCA_7125.1", "well formed GenBank"),
    ("refseq_id", "GCA_000007125.1", "well formed RefSeq"),
    ("assembly_name", "", "assembly name is empty"),
    ("strain", "tes;tus", "semicolon"),
])
def test_metadata_that_breaks_the_rules_is_refused(dataset, scheme, field, value, expected,
                                                   caplog):
    scores, outdir = scheme
    metadata = {**METADATA, field: value}
    assert export_main(scores, dataset.hierarchy, dataset.reference, outdir, **metadata) == 2
    assert expected in caplog.text


def test_a_score_file_from_another_tool_is_refused(dataset, scheme, tmp_path):
    scores, outdir = scheme
    wrong = tmp_path / "wrong.csv"
    wrong.write_text("position,group\n100,A.1\n")
    assert export_main(str(wrong), dataset.hierarchy, dataset.reference, outdir, **METADATA) == 1


def test_only_canSNPs_that_passed_are_exported(dataset, scheme, tmp_path, caplog):
    # A canSNP evaluate marked as failing is kept in its output but is not a
    # candidate for the scheme, so its node is left without one.
    scores, outdir = scheme
    failing = vcf_rewrite(scores, tmp_path / "failing.vcf",
                          change=lambda fields: fields[:6] + (
                              ["LowScore"] if fields[:2] == ["c1", "601"] else [fields[6]])
                          + fields[7:])
    assert export_main(failing, dataset.hierarchy, dataset.reference, outdir, **METADATA) == 2
    assert "1 nodes have no canSNP: B.1" in caplog.text


# -- against the real validator ----------------------------------------------


CANSNPER2_DATA = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "CanSNPer2-data", "validate.py")


@pytest.mark.skipif(not os.path.exists(CANSNPER2_DATA),
                    reason="CanSNPer2-data is not checked out beside this repository")
def test_the_exported_scheme_passes_the_cansnper2_data_validator(dataset, scheme, tmp_path):
    scores, outdir = scheme
    assert export_main(scores, dataset.hierarchy, dataset.reference, outdir, **METADATA) == 0

    # That script takes the repository layout from its own location, so give it
    # one holding nothing but the exported scheme.
    checkout = tmp_path / "checkout"
    (checkout / "source").mkdir(parents=True)
    shutil.copy(CANSNPER2_DATA, checkout / "validate.py")
    (checkout / "schemes.toml").write_text("")
    for name in os.listdir(outdir):
        shutil.copy(os.path.join(outdir, name), checkout / "source" / name)

    completed = subprocess.run([sys.executable, "validate.py"], cwd=str(checkout),
                               capture_output=True, text=True)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "All 1 schemes valid" in completed.stdout
