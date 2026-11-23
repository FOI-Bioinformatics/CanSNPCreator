"""The exit paths of the export subcommand's main function.

Each case pins the exit status and the records logged, so that main keeps
its messages and their order.
"""

import logging
import os

import pytest
from conftest import METADATA

from cansnpscorer.cansnp_export import main as export_main

NAMING = (
    "2 reference sequences are named neither by a versioned GenBank accession nor by a RefSeq "
    "accession of the NC_, NZ_ or AC_ kind, which would be replaced by its GenBank name; they are "
    "written as they are: c1, c2. CanSNPer2-download will not find a genome with sequences of these "
    "names, so the reference must be provided by hand"
)
NO_ROOT = (
    "the root 'root' has no canSNP, and a scheme cannot be written without one. A root's canSNP "
    "needs a genome outside the root in the alignment, an outgroup placed at no node (collapse "
    "--outgroup); with one, lower the evaluate cutoff, or root the scheme elsewhere"
)


def records(caplog: pytest.LogCaptureFixture) -> list[tuple[int, str]]:
    """The records of the export and of the report of its checks."""
    return [
        (record.levelno, record.getMessage())
        for record in caplog.records
        if record.module in ("cansnp_export", "validate")
    ]


def failing(scores: str, path: str, nodes: tuple[str, ...]) -> str:
    """A copy of the scores in which the canSNPs of these nodes failed."""
    with open(scores) as handle:
        lines = handle.read().splitlines()
    with open(path, "w") as handle:
        for line in lines:
            fields = line.split("\t")
            if not line.startswith("#") and any(f"NODE={node};" in fields[7] for node in nodes):
                fields[6] = "LowScore"
            handle.write("\t".join(fields) + "\n")
    return path


@pytest.mark.parametrize("given", ["genome_map", "alignment"])
def test_a_genome_map_and_an_alignment_go_together(written, dataset, tmp_path, caplog, given):
    _outdir, scores = written
    outdir = str(tmp_path / "out")
    options = {given: getattr(dataset, given)}
    caplog.set_level(logging.INFO)
    assert export_main(scores, dataset.hierarchy, dataset.reference, outdir, **METADATA, **options) == 1
    assert [(record.levelno, record.getMessage()) for record in caplog.records] == [
        (
            logging.ERROR,
            "the genomes file needs both --genome-map and --alignment: where each genome sits, "
            "and which genomes were aligned",
        )
    ]
    assert not os.path.exists(outdir)


def test_an_unreadable_score_file_exits_one(dataset, tmp_path, caplog):
    outdir = str(tmp_path / "out")
    status = export_main(
        str(tmp_path / "missing.vcf"), dataset.hierarchy, dataset.reference, outdir, **METADATA
    )
    assert status == 1
    found = [record.getMessage() for record in caplog.records if record.levelno >= logging.ERROR]
    assert len(found) == 1 and "missing.vcf" in found[0]
    assert not os.path.exists(outdir)


def test_nodes_without_a_canSNP_stop_the_export(written, dataset, tmp_path, caplog):
    _outdir, scores = written
    high = failing(scores, str(tmp_path / "high.vcf"), ("root", "A.1", "B"))
    outdir = str(tmp_path / "out")
    caplog.set_level(logging.INFO)
    assert export_main(high, dataset.hierarchy, dataset.reference, outdir, **METADATA) == 2
    assert records(caplog) == [
        (logging.WARNING, NAMING),
        (logging.ERROR, NO_ROOT),
        (
            logging.ERROR,
            "3 nodes have no canSNP: root, A.1, B. Every node of a CanSNPer2 scheme carries exactly "
            "one, so either lower the evaluate cutoff or pass --prune to leave these nodes out",
        ),
        (logging.ERROR, "2 problems with the scheme; nothing was written"),
    ]
    assert not os.path.exists(outdir)


def test_pruning_does_not_take_away_the_root(written, dataset, tmp_path, caplog):
    _outdir, scores = written
    high = failing(scores, str(tmp_path / "high.vcf"), ("root", "A.1", "B"))
    outdir = str(tmp_path / "out")
    caplog.set_level(logging.INFO)
    status = export_main(
        high,
        dataset.hierarchy,
        dataset.reference,
        outdir,
        **METADATA,
        prune=True,
        genome_map=dataset.genome_map,
        alignment=dataset.alignment,
    )
    assert status == 2
    assert records(caplog) == [
        (logging.WARNING, NAMING),
        (
            logging.WARNING,
            "3 nodes have no canSNP and were left out of the scheme, so their children now hang from "
            "the nearest node that has one: root, A.1, B. 2 genomes of the genome map sat at these "
            "nodes and now sit at an ancestor; the deepest pruned depth is 3 (the root is depth 1)",
        ),
        (logging.ERROR, NO_ROOT),
        (logging.ERROR, "1 problems with the scheme; nothing was written"),
    ]


def test_a_pruned_export_with_a_canSNP_set_and_genomes(written, dataset, tmp_path, caplog):
    _outdir, scores = written
    high = failing(scores, str(tmp_path / "high.vcf"), ("A.1", "B"))
    outdir = str(tmp_path / "out")
    caplog.set_level(logging.INFO)
    status = export_main(
        high,
        dataset.hierarchy,
        dataset.reference,
        outdir,
        **METADATA,
        prune=True,
        genome_map=dataset.genome_map,
        alignment=dataset.alignment,
        snps_per_node=2,
    )
    assert status == 0
    files = [
        os.path.join(outdir, "testus_organismus" + suffix)
        for suffix in ("_tree.txt", "_references.txt", "_canSNPs.txt", "_genomes.txt")
    ]
    assert records(caplog) == [
        (logging.WARNING, NAMING),
        (
            logging.WARNING,
            "2 nodes have no canSNP and were left out of the scheme, so their children now hang from "
            "the nearest node that has one: A.1, B. 2 genomes of the genome map sat at these nodes "
            "and now sit at an ancestor; the deepest pruned depth is 3 (the root is depth 1)",
        ),
        (logging.INFO, "No problems found in the scheme, with 2 warnings"),
        (logging.INFO, "The canSNP set holds 5 canSNPs, 1 to 1 per node"),
        (
            logging.WARNING,
            "5 nodes have fewer than 2 canSNPs that passed and lie 1000 bases apart: root, A, A.2, B.1, B.2",
        ),
        (
            logging.WARNING,
            "2 genomes sit at an inner node and were left out of the genomes file, which lists the "
            "genomes of the terminal nodes; they were used to find the canSNPs: s1 (A), s1b (A)",
        ),
        (
            logging.WARNING,
            "4 genomes hold no assembly accession in their names and were left out of the genomes "
            "file: ref, s2, s3, s4",
        ),
        (logging.INFO, "The genomes file lists 0 representatives, 0 members, 0 outgroups"),
        (logging.INFO, f"Wrote a scheme of 5 nodes: {', '.join(files)}"),
        (
            logging.INFO,
            "Check it with the CanSNPer2-data repository's validate.py before building a database",
        ),
    ]


def test_the_files_of_an_earlier_export_are_removed(written, dataset, tmp_path, caplog):
    _outdir, scores = written
    outdir = tmp_path / "out"
    assert (
        export_main(
            scores,
            dataset.hierarchy,
            dataset.reference,
            str(outdir),
            **METADATA,
            genome_map=dataset.genome_map,
            alignment=dataset.alignment,
        )
        == 0
    )
    (outdir / "testus_organismus_snp.txt").write_text("old\n")
    caplog.clear()
    caplog.set_level(logging.INFO)
    assert export_main(scores, dataset.hierarchy, dataset.reference, str(outdir), **METADATA) == 0
    removed = [message for level, message in records(caplog) if message.startswith("Removed")]
    assert removed == [
        f"Removed {outdir / name}, which an earlier export wrote and this one does not"
        for name in ("testus_organismus_snp.txt", "testus_organismus_genomes.txt")
    ]
    assert sorted(os.listdir(outdir)) == [
        "testus_organismus_canSNPs.txt",
        "testus_organismus_references.txt",
        "testus_organismus_tree.txt",
    ]
