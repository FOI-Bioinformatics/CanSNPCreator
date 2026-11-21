"""Export of a scored scheme as CanSNPer2-data source files.

The rules checked here are the ones that repository's schemes/README.md sets
out and its validate.py enforces. Some tests run that validate.py itself, from
the checkout CANSNPER2_DATA names or one beside this repository, and one builds
the exported scheme with CanSNPer2 2.3.0, from the python CANSNPER2_PYTHON
names or this one.
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
from cansnpscorer.helpers.scoring_methods import DEFAULT_WEIGHTS, TERMS, ScoredSnp

METADATA = {
    "organism": "testus_organismus",
    "genome": "TEST1",
    "strain": "testus",
    "citation": "Sjodin2026",
    "genbank_id": "GCA_000007125.1",
    "refseq_id": "GCF_000007125.1",
    "assembly_name": "ASM712v1",
}


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
                   "--genes", dataset.gff, "--repeats", dataset.repeats,
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


def test_positions_are_given_per_sequence(dataset, scheme):
    # From CanSNPer2 2.3.0 each SNP names its reference sequence and gives a
    # 1-based position on it, so a SNP on the second sequence is not offset
    # by the length of the first.
    scores, outdir = scheme
    assert export_main(scores, dataset.hierarchy, dataset.reference, outdir, **METADATA) == 0
    lines = read(os.path.join(outdir, "testus_organismus_snp.txt")).splitlines()
    assert lines[0].split("\t") == ["snp_id", "strain", "reference", "genome", "sequence",
                                    "sequence_position", "derived_base", "ancestral_base"]
    positions = {}
    for line in lines[1:]:
        fields = line.split("\t")
        positions[fields[0]] = (fields[4], int(fields[5]))

    for node, (column, _ancestral, _derived) in PLANTED.items():
        assert positions[node] == dataset.contig_of(column), node
    assert positions["B.2"] == ("c2", 1), "the first base of the second sequence"


def test_the_polarity_survives_the_round_trip(dataset, scheme):
    scores, outdir = scheme
    assert export_main(scores, dataset.hierarchy, dataset.reference, outdir, **METADATA) == 0
    for line in read(os.path.join(outdir, "testus_organismus_snp.txt")).splitlines()[1:]:
        node, _strain, _reference, _genome, _sequence, _position, derived, ancestral = (
            line.split("\t"))
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


def test_prune_reports_genomes_moved(dataset, scheme, tmp_path, caplog):
    scores, outdir = scheme
    trimmed = vcf_rewrite(scores, tmp_path / "fewer.vcf",
                          keep=lambda fields: fields[:2] != ["c1", "201"])
    genome_map = tmp_path / "pruned_map.tsv"
    genome_map.write_text("ref\tA\ns1\tA\ns3\tB.1\n")

    assert export_main(trimmed, dataset.hierarchy, dataset.reference, outdir, prune=True,
                       genome_map=str(genome_map), alignment=dataset.alignment, **METADATA) == 0
    assert "2 genomes of the genome map sat at these nodes and now sit at an ancestor" \
        in caplog.text
    assert "deepest pruned depth is 2" in caplog.text


def test_prune_without_a_genome_map_keeps_its_warning(dataset, scheme, tmp_path, caplog):
    scores, outdir = scheme
    trimmed = vcf_rewrite(scores, tmp_path / "fewer.vcf",
                          keep=lambda fields: fields[:2] != ["c1", "201"])

    assert export_main(trimmed, dataset.hierarchy, dataset.reference, outdir,
                       prune=True, **METADATA) == 0
    assert "genome map sat" not in caplog.text


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
    _scores, outdir = scheme
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
#
# CanSNPer2-data's validate.py takes the repository layout from its own
# location, so it is copied into a checkout holding nothing but the exported
# scheme, in schemes/<organism>/ with a README in the required layout. The
# checkout is taken from CANSNPER2_DATA, or found beside this repository.


CANSNPER2_DATA = os.environ.get("CANSNPER2_DATA") or os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "CanSNPer2-data")
VALIDATOR = os.path.join(CANSNPER2_DATA, "validate.py")


def validator_skip_reason():
    if not os.path.exists(VALIDATOR):
        return f"no CanSNPer2-data checkout at {CANSNPER2_DATA}; set CANSNPER2_DATA"
    if '"sequence_position"' not in read(VALIDATOR):
        return (f"the validate.py of {CANSNPER2_DATA} expects the SNP file of CanSNPer2 "
                f"before 2.3.0; set CANSNPER2_DATA to a checkout with the 2.3.0 format")
    return None


needs_validator = pytest.mark.skipif(validator_skip_reason() is not None,
                                     reason=validator_skip_reason() or "")


def scheme_checkout(outdir, tmp_path, organism="testus_organismus"):
    """A copy of the validator with the exported scheme as its only scheme."""
    checkout = tmp_path / "checkout"
    folder = checkout / "schemes" / organism
    folder.mkdir(parents=True)
    shutil.copy(VALIDATOR, checkout / "validate.py")
    (checkout / "schemes.toml").write_text("")
    for name in os.listdir(outdir):
        shutil.copy(os.path.join(outdir, name), folder / name)

    nodes = len(read(os.path.join(outdir, f"{organism}_tree.txt")).splitlines())
    per_node = "1"
    snp_set_path = os.path.join(outdir, f"{organism}_canSNPs.txt")
    if os.path.exists(snp_set_path):
        ranks = [int(line.split("\t")[1]) for line in read(snp_set_path).splitlines()[1:]]
        if max(ranks) > 1:
            per_node = f"up to {max(ranks)}"
    facts = [("Scheme", f"`{organism}`"), ("Status", "draft"), ("Organisms", "Testus"),
             ("Reference genome", METADATA["genbank_id"]), ("Nodes", str(nodes)),
             ("canSNPs per node", per_node), ("Origin", "CanSNPCreator test"),
             ("Last changed", "2026-11-21")]
    sections = ["Description", "Origin", "Validation", "Known issues", "References",
                "Changes"]
    (folder / "README.md").write_text(
        f"# {organism}\n\n| | |\n|---|---|\n"
        + "".join(f"| {key} | {value} |\n" for key, value in facts)
        + "".join(f"\n## {section}\n\nNone.\n" for section in sections))
    return checkout


def run_validator(checkout):
    return subprocess.run([sys.executable, "validate.py"], cwd=str(checkout),
                          capture_output=True, text=True)


@needs_validator
def test_the_exported_scheme_passes_the_cansnper2_data_validator(dataset, scheme, tmp_path):
    scores, outdir = scheme
    assert export_main(scores, dataset.hierarchy, dataset.reference, outdir, **METADATA) == 0

    completed = run_validator(scheme_checkout(outdir, tmp_path))
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "All 1 schemes valid" in completed.stdout


# -- names of the reference sequences ------------------------------------------
#
# CanSNPer2 2.3.0 names each reference sequence by its GenBank accession, as
# in the FASTA files CanSNPer2-download fetches. A reference whose sequences
# carry RefSeq names is written with their GenBank names from an alias table.

REFSEQ_NAMES = {"c1": "NC_000001.1", "c2": "NC_000002.1"}
GENBANK_NAMES = {"c1": "CP000001.1", "c2": "CP000002.1"}


@pytest.fixture
def refseq(dataset, scheme, tmp_path):
    """The scores and the reference of the dataset, with RefSeq sequence names,
    and an NCBI assembly report pairing them with GenBank names."""
    from dataclasses import replace

    from conftest import write_fasta

    from cansnpscorer import vcf
    from cansnpscorer.helpers.file_utils import reference_parser
    from cansnpscorer.model import Reference

    scores, outdir = scheme
    reference = reference_parser(dataset.reference)
    renamed = Reference.from_records(
        (REFSEQ_NAMES[contig.name], contig.sequence) for contig in reference)
    fasta = str(tmp_path / "refseq.fasta")
    write_fasta(fasta, [(contig.name, contig.sequence) for contig in renamed])
    planted = [replace(score, contig=REFSEQ_NAMES[score.contig])
               for score in vcf.read_scores(scores, reference)]
    path = str(tmp_path / "refseq.vcf")
    vcf.write_scores(path, planted, renamed, fasta, 0, None, DEFAULT_WEIGHTS)

    report = tmp_path / "assembly_report.txt"
    report.write_text(
        "# Assembly name:  ASM712v1\n"
        "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type"
        "\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length"
        "\tUCSC-style-name\n"
        + "".join(f"chr{index}\tassembled-molecule\t{index}\tChromosome\t"
                  f"{GENBANK_NAMES[name]}\t=\t{REFSEQ_NAMES[name]}\tPrimary Assembly\t"
                  f"{reference.contig(name).length}\tna\n"
                  for index, name in enumerate(("c1", "c2"), 1)))
    return path, fasta, str(report), outdir


def sequences_of(outdir):
    lines = read(os.path.join(outdir, "testus_organismus_snp.txt")).splitlines()[1:]
    return {fields[0]: (fields[4], int(fields[5]))
            for fields in (line.split("\t") for line in lines)}


def test_refseq_contigs_are_written_with_genbank_names(dataset, refseq):
    scores, fasta, report, outdir = refseq
    metadata = [f"--{key.replace('_', '-')}={value}" for key, value in METADATA.items()]
    assert cli_run("export", "--scores", scores, "--hierarchy", dataset.hierarchy,
                   "--reference", fasta, "--outdir", outdir, *metadata,
                   "--contig-aliases", report, "--snps-per-node", "2") == 0

    found = sequences_of(outdir)
    for node, (column, _ancestral, _derived) in PLANTED.items():
        contig, position = dataset.contig_of(column)
        assert found[node] == (GENBANK_NAMES[contig], position), node
    snp_set = read(os.path.join(outdir, "testus_organismus_canSNPs.txt"))
    assert "NC_" not in snp_set
    assert {line.split("\t")[3] for line in snp_set.splitlines()[1:]} \
        == set(GENBANK_NAMES.values())


def test_refseq_contigs_without_aliases_are_refused(dataset, refseq, caplog):
    scores, fasta, _report, outdir = refseq
    assert export_main(scores, dataset.hierarchy, fasta, outdir, **METADATA) == 2
    assert "NC_000001.1" in caplog.text
    assert "--contig-aliases <assembly_report.txt>" in caplog.text
    assert not os.path.exists(os.path.join(outdir, "testus_organismus_snp.txt"))


def test_refseq_contigs_with_no_genbank_alias_are_refused(dataset, refseq, tmp_path, caplog):
    scores, fasta, _report, outdir = refseq
    aliases = tmp_path / "chromAlias.txt"
    aliases.write_text("NC_000001.1\tchr1\nNC_000002.1\tchr2\n")
    assert export_main(scores, dataset.hierarchy, fasta, outdir, contig_aliases=str(aliases),
                       **METADATA) == 2
    assert "NC_000002.1" in caplog.text
    # The file was given; asking for --contig-aliases again would not help.
    refusals = [record.getMessage() for record in caplog.records
                if "NC_000002.1" in record.getMessage()]
    assert any(f"{aliases} gives no GenBank accession" in message for message in refusals)
    assert not any("give it with --contig-aliases" in message for message in refusals)


def test_a_refseq_sequence_not_identical_to_its_genbank_one_is_refused(dataset, refseq,
                                                                       caplog):
    # NCBI gives the Relationship <> when the RefSeq sequence is not identical
    # to the GenBank one, so a position on one need not hold on the other.
    scores, fasta, report, outdir = refseq
    with open(report) as handle:
        text = handle.read()
    with open(report, "w") as handle:
        handle.write(text.replace("CP000002.1\t=\t", "CP000002.1\t<>\t"))
    assert export_main(scores, dataset.hierarchy, fasta, outdir, contig_aliases=report,
                       **METADATA) == 2
    refusals = [record.getMessage() for record in caplog.records
                if "not identical" in record.getMessage()]
    assert len(refusals) == 1
    assert "NC_000002.1" in refusals[0] and "CP000002.1" in refusals[0]
    assert "NC_000001.1" not in refusals[0]
    assert not os.path.exists(os.path.join(outdir, "testus_organismus_snp.txt"))


def test_the_alias_reader_keeps_which_sequences_are_not_identical(tmp_path):
    from cansnpscorer.helpers.file_utils import read_contig_aliases

    report = tmp_path / "assembly_report.txt"
    report.write_text(
        "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type"
        "\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length"
        "\tUCSC-style-name\n"
        "chr1\tassembled-molecule\t1\tChromosome\tCP000001.1\t=\tNC_000001.1\tPrimary"
        "\t10\tna\n"
        "pX\tassembled-molecule\tpX\tPlasmid\tCP000002.1\t<>\tNC_000002.1\tPrimary"
        "\t10\tna\n")
    same, different = read_contig_aliases(str(report))
    assert same == {"chr1", "CP000001.1", "NC_000001.1"} and not same.differs
    assert different == {"pX", "CP000002.1", "NC_000002.1"}
    assert different.differs == {"NC_000002.1"}


@pytest.mark.parametrize("name", ["NT_000003.1", "NW_000004.1"])
def test_a_refseq_contig_or_scaffold_name_is_not_called_neither_kind(name):
    from cansnpscorer.cansnp_export import sequence_names

    names, problems, warnings = sequence_names([name, "contig_7"], None)
    assert names == {name: name, "contig_7": "contig_7"} and problems == []
    (warning,) = warnings
    assert "neither a GenBank nor a RefSeq accession" not in warning
    assert name in warning and "contig_7" in warning
    assert "NT_" in warning and "NW_" in warning


def test_a_local_sequence_name_is_written_with_one_warning(dataset, scheme, caplog):
    scores, outdir = scheme
    assert export_main(scores, dataset.hierarchy, dataset.reference, outdir, **METADATA) == 0
    warnings = [record.getMessage() for record in caplog.records
                if record.levelname == "WARNING" and "CanSNPer2-download" in record.getMessage()]
    assert len(warnings) == 1
    assert "c1" in warnings[0] and "c2" in warnings[0]
    assert sequences_of(outdir)["B.2"] == ("c2", 1)


def test_an_unreadable_alias_file_is_refused(dataset, refseq, tmp_path, caplog):
    scores, fasta, _report, outdir = refseq
    empty = tmp_path / "empty.txt"
    empty.write_text("")
    assert export_main(scores, dataset.hierarchy, fasta, outdir, contig_aliases=str(empty),
                       **METADATA) == 1
    assert "empty.txt" in caplog.text


# -- against CanSNPer2 2.3.0 ---------------------------------------------------


CANSNPER2_PYTHON = os.environ.get("CANSNPER2_PYTHON") or sys.executable


def cansnper2_230_missing():
    try:
        completed = subprocess.run(
            [CANSNPER2_PYTHON, "-c", "import CanSNPer2.modules.SchemeInfo"],
            capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return True
    return completed.returncode != 0


@pytest.mark.skipif(cansnper2_230_missing(),
                    reason="CanSNPer2 2.3.0 is not importable; set CANSNPER2_PYTHON to a python "
                           "that has it")
def test_the_export_builds_with_cansnper2_230(dataset, scheme, tmp_path):
    scores, outdir = scheme
    assert export_main(scores, dataset.hierarchy, dataset.reference, outdir, **METADATA,
                       snps_per_node=2, spacing=0) == 0
    database = str(tmp_path / "testus.db")
    completed = subprocess.run(
        [CANSNPER2_PYTHON, "-m", "CanSNPer2.SNPDatabase", "--create", "-db", database,
         "--tree", os.path.join(outdir, "testus_organismus_tree.txt"),
         "--annotation", os.path.join(outdir, "testus_organismus_snp.txt"),
         "--references", os.path.join(outdir, "testus_organismus_references.txt"),
         "--tmpdir", str(tmp_path)],
        cwd=str(tmp_path), capture_output=True, text=True, timeout=300)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert os.path.getsize(database) > 0

    # The database gives back the SNPs as they were written.
    exported = tmp_path / "exported"
    completed = subprocess.run(
        [CANSNPER2_PYTHON, "-m", "CanSNPer2.SNPDatabase", "-db", database, "--export",
         "-o", str(exported)],
        cwd=str(tmp_path), capture_output=True, text=True, timeout=300)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    written = sequences_of(outdir)
    found = {}
    for path in exported.rglob("*_snp*.txt"):
        rows = [line.split("\t") for line in path.read_text().splitlines()]
        at = {name: index for index, name in enumerate(rows[0])}
        for row in rows[1:]:
            found[row[at["snp_id"]]] = (row[at["sequence"]], int(row[at["sequence_position"]]))
    assert found == written


# -- more than one canSNP per node -------------------------------------------


@pytest.fixture
def several(dataset, scheme, tmp_path):
    """A score file with, for each planted node, two more canSNPs: one 5 bases
    from the planted one and scoring just below it, one 20 bases away and
    scoring lower still."""
    from cansnpscorer import vcf
    from cansnpscorer.helpers.file_utils import reference_parser

    scores, outdir = scheme
    reference = reference_parser(dataset.reference)
    planted = vcf.read_scores(scores, reference)
    extra = []
    for original in planted:
        contig = reference.contig(original.contig)
        for offset, drop in ((5, 1), (20, 2)):
            position = original.position + offset
            base = contig.base(position - 1)
            derived = "A" if base != "A" else "C"
            extra.append(ScoredSnp(contig=original.contig, position=position, node=original.node,
                                   total=original.total - drop, terms=dict(original.terms),
                                   derived=derived, ancestral=base))
    path = str(tmp_path / "several.vcf")
    vcf.write_scores(path, planted + extra, reference, dataset.reference, 0, None,
                     DEFAULT_WEIGHTS)
    return path, outdir


def snp_set(outdir):
    rows = [line.split("\t") for line in
            read(os.path.join(outdir, "testus_organismus_canSNPs.txt")).splitlines()]
    assert rows[0] == ["node", "rank", "genome", "sequence", "sequence_position",
                       "derived_base", "ancestral_base", "score"]
    return rows[1:]


def test_one_snp_per_node_writes_no_snp_set(dataset, several):
    scores, outdir = several
    assert export_main(scores, dataset.hierarchy, dataset.reference, outdir, **METADATA) == 0
    assert not os.path.exists(os.path.join(outdir, "testus_organismus_canSNPs.txt"))


def test_the_snp_set_ranks_each_node_and_keeps_the_snp_file_as_rank_one(dataset, several):
    scores, outdir = several
    assert export_main(scores, dataset.hierarchy, dataset.reference, outdir, **METADATA,
                       snps_per_node=None, spacing=0) == 0
    rows = snp_set(outdir)
    assert len(rows) == 3 * len(PLANTED)
    snp = [line.split("\t") for line in
           read(os.path.join(outdir, "testus_organismus_snp.txt")).splitlines()[1:]]
    first = {row[0]: (row[3], row[4], row[5], row[6]) for row in rows if row[1] == "1"}
    assert first == {row[0]: (row[4], row[5], row[6], row[7]) for row in snp}
    for node in PLANTED:
        ranked = [row for row in rows if row[0] == node]
        assert [row[1] for row in ranked] == ["1", "2", "3"]
        assert [int(row[7]) for row in ranked] == sorted((int(row[7]) for row in ranked),
                                                         reverse=True)


def test_spacing_passes_over_a_canSNP_beside_one_already_kept(dataset, several):
    scores, outdir = several
    assert export_main(scores, dataset.hierarchy, dataset.reference, outdir, **METADATA,
                       snps_per_node=2, spacing=10) == 0
    rows = snp_set(outdir)
    for node in PLANTED:
        positions = [int(row[4]) for row in rows if row[0] == node]
        # The one 5 bases away gives way to the one 20 bases away.
        assert len(positions) == 2 and positions[1] - positions[0] == 20


def test_the_snp_set_is_chosen_from_the_command_line(dataset, several):
    scores, outdir = several
    metadata = [f"--{key.replace('_', '-')}={value}" for key, value in METADATA.items()]
    assert cli_run("export", "--scores", scores, "--hierarchy", dataset.hierarchy,
                   "--reference", dataset.reference, "--outdir", outdir, *metadata,
                   "--snps-per-node", "all", "--min-spacing", "0") == 0
    assert len(snp_set(outdir)) == 3 * len(PLANTED)


@needs_validator
def test_the_snp_set_passes_the_cansnper2_data_validator(dataset, several, tmp_path):
    scores, outdir = several
    assert export_main(scores, dataset.hierarchy, dataset.reference, outdir, **METADATA,
                       snps_per_node=None, spacing=0) == 0
    checkout = scheme_checkout(outdir, tmp_path)
    completed = run_validator(checkout)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "All 1 schemes valid" in completed.stdout

    # A rank 1 that differs from the SNP file is refused, so the file is read.
    snp_set_path = checkout / "schemes" / "testus_organismus" / "testus_organismus_canSNPs.txt"
    lines = snp_set_path.read_text().splitlines()
    lines[1], lines[2] = lines[2].replace("\t2\t", "\t1\t", 1), lines[1].replace("\t1\t", "\t2\t", 1)
    snp_set_path.write_text("\n".join(lines) + "\n")
    completed = run_validator(checkout)
    assert completed.returncode == 1
    assert "is not the SNP" in completed.stdout


# -- the genomes the scheme was made from ------------------------------------


def test_genome_rows_give_each_genome_its_accession_node_and_role():
    from cansnpscorer.cansnp_export import genome_rows
    from cansnpscorer.hierarchy import Hierarchy

    tree = Hierarchy.from_child_parent([("A", "root"), ("A.1", "A"), ("B", "root")])
    clusters = {
        "A.1": ["Fam_Gen_sp_GCF_000000001.1", "GCA_000000002.1"],
        "B": ["Fam_Gen_sp_GCF_000000003.1", "in-house-7"],
    }
    aligned = ["Fam_Gen_sp_GCF_000000001.1", "Fam_Gen_sp_GCF_000000003.1",
               "Fam_Out_grp_GCF_000000009.1"]
    # A.1 was left out of the scheme, so its genomes move to A.
    rows, warnings = genome_rows(tree, ["root", "A", "B"], clusters, aligned)
    assert rows == [
        ("GCF_000000001.1", "A", "representative"),
        ("GCF_000000003.1", "B", "representative"),
        ("GCA_000000002.1", "A", "member"),
        ("GCF_000000009.1", "", "outgroup"),
    ]
    assert "in-house-7" in warnings[0]


def accessioned(dataset, tmp_path):
    """The planted dataset with RepGenR2-style genome names."""
    from conftest import CLUSTERS, GENOMES, write_fasta

    from cansnpscorer.helpers.file_utils import read_alignment

    names = {genome: f"Fam_Gen_sp_GCF_00000000{index}.1" for index, genome
             in enumerate(GENOMES, 1)}
    alignment = read_alignment(dataset.alignment)
    aligned = str(tmp_path / "aligned.fasta")
    write_fasta(aligned, [(names[name], sequence) for name, sequence in alignment.items()])
    genome_map = tmp_path / "genome_map.tsv"
    genome_map.write_text("".join(f"{names[genome]}\t{node}\n" for node, genomes
                                  in CLUSTERS.items() for genome in genomes)
                          + "GCA_000000099.1\tA.1\n")  # a member, not aligned
    return aligned, str(genome_map)


def test_the_genomes_file_is_written_validated_and_audited(dataset, tmp_path):
    aligned, genome_map = accessioned(dataset, tmp_path)
    snps, scores = str(tmp_path / "snps.vcf"), str(tmp_path / "scores.vcf")
    outdir = str(tmp_path / "scheme")
    assert cli_run("extract", "--alignment", aligned, "--reference", dataset.reference,
                   "--hierarchy", dataset.hierarchy, "--genome-map", genome_map,
                   "-o", snps) == 0
    assert cli_run("evaluate", "--snps", snps, "--reference", dataset.reference,
                   "--genes", dataset.gff, "--repeats", dataset.repeats,
                   "-o", scores, "--cutoff", "0") == 0
    assert export_main(scores, dataset.hierarchy, dataset.reference, outdir,
                       genome_map=genome_map, alignment=aligned, **METADATA) == 0

    lines = read(os.path.join(outdir, "testus_organismus_genomes.txt")).splitlines()
    assert lines[0] == "genome\tnode\trole"
    rows = [line.split("\t") for line in lines[1:]]
    assert ["GCF_000000002.1", "A.1", "representative"] in rows
    assert ["GCA_000000099.1", "A.1", "member"] in rows
    assert ["GCF_000000001.1", "", "outgroup"] in rows, "the reference is placed at no node"

    # audit uses the scheme's own genomes file when given only the alignment
    report = tmp_path / "findings.tsv"
    from cansnpscorer.cansnp_audit import main as audit_main
    assert audit_main(outdir, reference=dataset.reference, alignment=aligned,
                      output=str(report)) == 0
    assert "not-canonical" not in report.read_text()
    assert "single-site" in report.read_text(), "the genome checks ran"


@needs_validator
def test_the_genomes_file_passes_the_cansnper2_data_validator(dataset, tmp_path):
    aligned, genome_map = accessioned(dataset, tmp_path)
    snps, scores = str(tmp_path / "snps.vcf"), str(tmp_path / "scores.vcf")
    outdir = str(tmp_path / "scheme")
    assert cli_run("extract", "--alignment", aligned, "--reference", dataset.reference,
                   "--hierarchy", dataset.hierarchy, "--genome-map", genome_map,
                   "-o", snps) == 0
    assert cli_run("evaluate", "--snps", snps, "--reference", dataset.reference,
                   "--genes", dataset.gff, "--repeats", dataset.repeats,
                   "-o", scores, "--cutoff", "0") == 0
    assert export_main(scores, dataset.hierarchy, dataset.reference, outdir,
                       genome_map=genome_map, alignment=aligned, **METADATA) == 0

    completed = run_validator(scheme_checkout(outdir, tmp_path))
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "All 1 schemes valid" in completed.stdout


def test_the_region_is_set_from_the_command_line(tmp_path, dataset, several):
    scores, outdir = several
    assert export_main(scores, dataset.hierarchy, dataset.reference, outdir, **METADATA,
                       snps_per_node=3, spacing=0) == 0
    found = tmp_path / "found.tsv"
    cli_run("audit", "--scheme-dir", outdir, "--set-region", "1000000", "-o", str(found))
    assert "clustered-set" in found.read_text()
    cli_run("audit", "--scheme-dir", outdir, "--set-region", "0", "-o", str(found))
    assert "clustered-set" not in found.read_text()
