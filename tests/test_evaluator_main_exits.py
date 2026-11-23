"""The exit paths of the evaluate subcommand's main function.

Each case pins the exit status and the records logged, so that main keeps
its messages and their order.
"""

import logging
import os

import pytest
from conftest import cli_run

from cansnpscorer import cansnp_evaluator
from cansnpscorer.cansnp_evaluator import main as evaluate_main


def own(caplog: pytest.LogCaptureFixture) -> list[tuple[int, str]]:
    """The records the evaluator module logs itself."""
    return [
        (record.levelno, record.getMessage())
        for record in caplog.records
        if record.module == "cansnp_evaluator"
    ]


def errors(caplog: pytest.LogCaptureFixture) -> list[str]:
    return [record.getMessage() for record in caplog.records if record.levelno >= logging.ERROR]


@pytest.fixture
def evaluate_inputs(dataset, tmp_path):
    """The arguments of a run of evaluate on the extracted SNPs of the dataset."""
    snp_file = str(tmp_path / "snps.vcf")
    assert (
        cli_run(
            "extract",
            "--genome-map",
            dataset.genome_map,
            "--hierarchy",
            dataset.hierarchy,
            "--alignment",
            dataset.alignment,
            "--reference",
            dataset.reference,
            "-o",
            snp_file,
        )
        == 0
    )
    return {
        "snps": snp_file,
        "reference": dataset.reference,
        "genes": dataset.gff,
        "repeats": dataset.repeats,
        "output": str(tmp_path / "scores.vcf"),
    }


def test_an_unreadable_snp_file_exits_one(evaluate_inputs, tmp_path, caplog):
    missing = str(tmp_path / "missing.vcf")
    assert evaluate_main(**{**evaluate_inputs, "snps": missing}) == 1
    assert len(errors(caplog)) == 1 and "missing.vcf" in errors(caplog)[0]
    assert not os.path.exists(evaluate_inputs["output"])


def test_a_gene_file_on_other_sequences_exits_one(evaluate_inputs, tmp_path, caplog):
    genes = tmp_path / "other.gff"
    genes.write_text("##gff-version 3\nzz\tR\tgene\t1\t10\t.\t+\t.\tID=g\n")
    assert evaluate_main(**{**evaluate_inputs, "genes": str(genes)}) == 1
    assert len(errors(caplog)) == 1
    assert errors(caplog)[0].startswith(
        "none of the 1 sequences of the gene file (zz) is a contig of the reference"
    )


def test_an_annotation_outside_its_contig_exits_two(evaluate_inputs, tmp_path, caplog):
    repeats = tmp_path / "long.bed"
    repeats.write_text("c1\t0\t100000\n")
    assert evaluate_main(**{**evaluate_inputs, "repeats": str(repeats)}) == 2
    assert errors(caplog) == [
        "1 repeat intervals on 'c1' lie outside the contig, which is 700 bases long; the first is 1..100000",
        "1 problems with the annotations; nothing was written",
    ]
    assert not os.path.exists(evaluate_inputs["output"])


def test_snps_that_cannot_be_placed_exit_two(evaluate_inputs, monkeypatch, caplog):
    # read_snps refuses what locate_snps would, so the failure is made here.
    def refuse(_snps, _reference):
        raise ValueError("1 of 8 SNPs do not fit the reference genome: c1:1 is outside the contig")

    monkeypatch.setattr(cansnp_evaluator, "locate_snps", refuse)
    assert evaluate_main(**evaluate_inputs) == 2
    assert errors(caplog) == ["1 of 8 SNPs do not fit the reference genome: c1:1 is outside the contig"]
    assert not os.path.exists(evaluate_inputs["output"])


def test_the_warnings_of_a_run_come_in_order(evaluate_inputs, dataset, tmp_path, caplog):
    genes = tmp_path / "nocds.gff"
    with open(dataset.gff) as handle:
        genes.write_text("".join(line for line in handle if "\tCDS\t" not in line))
    masks = tmp_path / "masks.bed"
    masks.write_text("c1\t100\t300\n")
    caplog.set_level(logging.INFO)
    status = evaluate_main(**{**evaluate_inputs, "genes": str(genes), "masks": str(masks), "flank": 10000})
    assert status == 0
    output = evaluate_inputs["output"]
    assert own(caplog) == [
        (
            logging.WARNING,
            "--flank 10000 is larger than the 1000 bases up to which extract records the distance to "
            "the nearest variation (INFO/NEAREST), so variation further away is not known: a canSNP "
            "with none within 1000 bases gets 0.10 of the full neighbour term, and distances beyond "
            "1000 are not told apart",
        ),
        (
            logging.WARNING,
            "The GFF3 file holds no CDS features, so no reading frame is known and every SNP scores 0 "
            "for the mutation term",
        ),
        (
            logging.INFO,
            "0 canSNPs lie in sequence that occurs more than once in the reference (31-mers, both "
            "strands); within a node they are ranked after those in unique sequence",
        ),
        (logging.INFO, "2 canSNPs lie in masked regions"),
        (logging.INFO, "Scored 7 canSNPs out of 8 SNPs; 5 pass the cutoff of 0"),
        (logging.INFO, f"Wrote 7 scored canSNPs to {output}"),
    ]


def test_the_messages_on_an_observed_table_come_in_order(evaluate_inputs, written, dataset, tmp_path, caplog):
    from cansnpscorer.cansnp_audit import main as audit_main

    outdir, _scores = written
    prefix = str(tmp_path / "typed")
    audit_main(
        outdir,
        reference=dataset.reference,
        alignment=dataset.alignment,
        genome_map=dataset.genome_map,
        typing_report=prefix,
    )
    with open(prefix + "_snps.tsv") as handle:
        lines = handle.read().splitlines()
    header = lines[0].split("\t")
    rows = [line.split("\t") for line in lines[1:]]
    derived, ancestral = header.index("derived_base"), header.index("ancestral_base")
    for row in rows[:3]:
        row[derived], row[ancestral] = "N", "N"
    observed = tmp_path / "other_snps.tsv"
    observed.write_text("".join("\t".join(fields) + "\n" for fields in [header, *rows]))

    caplog.set_level(logging.INFO)
    caplog.clear()
    assert evaluate_main(**{**evaluate_inputs, "observed": str(observed), "copy_kmer": 0}) == 0
    assert [message for _level, message in own(caplog)] == [
        "4 of 7 canSNPs were observed in the typing report, 0 of them disagreeing with their node in "
        "some genome and 0 lacking a base in some; within a node these are ranked after the others",
        f"3 rows of {observed} name the node and position of a canSNP here with other alleles, and "
        f"were not used; the report may have been made against another reference genome",
        "3 canSNPs were not in the scheme the report typed with. Within a node they are ranked after "
        "the observed ones that held and before those that did not; for a ranking of every "
        "candidate, type with a scheme exported with --snps-per-node all",
        "Scored 7 canSNPs out of 8 SNPs; 7 pass the cutoff of 0",
        f"Wrote 7 scored canSNPs to {evaluate_inputs['output']}",
    ]
