"""The exit paths of the audit subcommand's main function.

Each case pins the exit status and the error records logged, so that the
option checks keep their messages and their order.
"""

import logging

import pytest
from conftest import METADATA

from cansnpscorer.cansnp_audit import main as audit_main

READ = object()


def errors(caplog: pytest.LogCaptureFixture) -> list[str]:
    return [record.getMessage() for record in caplog.records if record.levelno >= logging.ERROR]


def own(caplog: pytest.LogCaptureFixture) -> list[tuple[int, str]]:
    """The records the audit module logs itself."""
    return [
        (record.levelno, record.getMessage()) for record in caplog.records if record.module == "cansnp_audit"
    ]


CHECKING = (
    "checking the scheme against the genomes needs --alignment, and either "
    "--genome-map or --tree to say where the genomes sit in the scheme, or the "
    "scheme's own genomes file"
)

#: Options, with READ standing for a path the test fills in, and the one error
#: the audit gives for them before it reads anything.
REFUSALS = [
    ({"scores": READ}, "--scores needs --reference, the genome the scores were made against"),
    (
        {"scores": READ, "masks": READ, "alignment": READ},
        "--scores needs --reference, the genome the scores were made against",
    ),
    ({"masks": READ}, "--masks needs --reference, which names the sequences the masks lie on"),
    (
        {"masks": READ, "alignment": READ},
        "--masks needs --reference, which names the sequences the masks lie on",
    ),
    ({"alignment": READ}, "--alignment needs --reference, the genome the alignment is against"),
    (
        {"alignment": READ, "genome": "TEST1"},
        "--alignment needs --reference, the genome the alignment is against",
    ),
    (
        {"reference": READ, "rule": "clade"},
        "--rule clade needs --tree, --genome-map with the scheme's genomes, and --outgroup",
    ),
    (
        {"reference": READ, "rule": "clade", "tree": READ, "genome_map": READ, "members": READ},
        "--rule clade needs --tree, --genome-map with the scheme's genomes, and --outgroup",
    ),
    ({"reference": READ, "outgroups": ("x",)}, "--outgroup applies to --rule clade"),
    ({"reference": READ, "outgroups": ("x",), "members": READ}, "--outgroup applies to --rule clade"),
    (
        {
            "reference": READ,
            "rule": "clade",
            "tree": READ,
            "genome_map": READ,
            "outgroups": ("x",),
            "members": READ,
            "alignment": READ,
        },
        "--members applies to --rule jaccard",
    ),
    ({"reference": READ, "genome_map": READ}, CHECKING),
    ({"reference": READ, "tree": READ}, CHECKING),
    ({"reference": READ, "genome_map": READ, "members": READ}, CHECKING),
    ({"reference": READ, "genome_map": READ, "tree": READ, "alignment": READ}, CHECKING),
    ({"reference": READ, "members": READ}, "--members applies to the leaves of --tree"),
    ({"reference": READ, "no_new_clades": True}, "--no-new-clades applies to --alignment"),
    (
        {"reference": READ, "no_new_clades": True, "typing_report": READ},
        "--no-new-clades applies to --alignment",
    ),
    (
        {"reference": READ, "typing_report": READ},
        "--typing-report types the genomes of --alignment, so it needs one",
    ),
    (
        {"typing_report": READ, "genome": "TEST1", "masks_concatenated": True},
        "--typing-report types the genomes of --alignment, so it needs one",
    ),
    (
        {"genome": "TEST1"},
        "--genome names the scheme's reference genome that --reference holds, so it needs --reference",
    ),
    (
        {"genome": "TEST1", "masks_concatenated": True},
        "--genome names the scheme's reference genome that --reference holds, so it needs --reference",
    ),
    ({"masks_concatenated": True}, "--masks-concatenated applies to --masks"),
]


@pytest.mark.parametrize(("options", "message"), REFUSALS)
def test_an_option_combination_is_refused_before_anything_is_read(written, dataset, caplog, options, message):
    outdir, _scores = written
    paths = {
        "scores": _scores,
        "masks": dataset.repeats,
        "alignment": dataset.alignment,
        "reference": dataset.reference,
        "tree": dataset.hierarchy_newick,
        "genome_map": dataset.genome_map,
        "members": dataset.genome_map,
        "typing_report": dataset.directory + "/typed",
    }
    given = {name: paths[name] if value is READ else value for name, value in options.items()}
    caplog.set_level(logging.INFO)
    assert audit_main(outdir, **given) == 1
    assert [(record.levelno, record.getMessage()) for record in caplog.records] == [(logging.ERROR, message)]


def test_a_scheme_that_cannot_be_read_exits_one(tmp_path, caplog):
    assert audit_main(str(tmp_path / "missing")) == 1
    assert len(errors(caplog)) == 1


def test_a_genome_the_scheme_does_not_annotate_exits_one(written, dataset, caplog):
    outdir, _scores = written
    assert audit_main(outdir, reference=dataset.reference, genome="ZZ") == 1
    assert errors(caplog) == [f"the scheme annotates no SNP on genome 'ZZ'; it names {METADATA['genome']}"]


def test_an_unreadable_score_file_exits_one(written, dataset, tmp_path, caplog):
    outdir, _scores = written
    assert audit_main(outdir, reference=dataset.reference, scores=str(tmp_path / "missing.vcf")) == 1
    assert len(errors(caplog)) == 1 and "missing.vcf" in errors(caplog)[0]


def test_an_alignment_with_no_placement_and_no_genomes_file_exits_one(written, dataset, caplog):
    outdir, _scores = written
    caplog.set_level(logging.INFO)
    assert audit_main(outdir, reference=dataset.reference, alignment=dataset.alignment) == 1
    assert own(caplog) == [
        (logging.INFO, "Reading only the scheme's columns of the alignment"),
        (
            logging.ERROR,
            "the scheme lists no genomes (testus_organismus_genomes.txt), so --genome-map or "
            "--tree must say where the genomes sit",
        ),
    ]


def test_an_unreadable_tree_exits_one_after_the_alignment_is_read(written, dataset, tmp_path, caplog):
    outdir, _scores = written
    caplog.set_level(logging.INFO)
    status = audit_main(
        outdir, reference=dataset.reference, alignment=dataset.alignment, tree=str(tmp_path / "missing.nwk")
    )
    assert status == 1
    assert own(caplog)[0] == (logging.INFO, "Reading every column of the alignment")
    assert len(errors(caplog)) == 1 and "missing.nwk" in errors(caplog)[0]


def test_a_typing_report_without_a_placement_types_only(written, dataset, tmp_path, caplog):
    outdir, _scores = written
    caplog.set_level(logging.INFO)
    status = audit_main(
        outdir,
        reference=dataset.reference,
        alignment=dataset.alignment,
        typing_report=str(tmp_path / "typed"),
    )
    assert status == 0
    messages = [message for _level, message in own(caplog)]
    assert messages[:2] == [
        "Reading only the scheme's columns of the alignment",
        "No genome map, tree or genomes file says where the genomes sit, so they are typed with the "
        "scheme but its SNPs are not checked against a placement",
    ]


def test_no_new_clades_reads_only_the_scheme_columns(written, dataset, caplog):
    outdir, _scores = written
    caplog.set_level(logging.INFO)
    status = audit_main(
        outdir,
        reference=dataset.reference,
        alignment=dataset.alignment,
        genome_map=dataset.genome_map,
        no_new_clades=True,
    )
    assert status == 0
    assert own(caplog)[0] == (logging.INFO, "Reading only the scheme's columns of the alignment")


def test_findings_that_cannot_be_written_exit_one(written, dataset, tmp_path, caplog):
    outdir, _scores = written
    assert audit_main(outdir, reference=dataset.reference, output=str(tmp_path / "no" / "f.tsv")) == 1
    assert "f.tsv" in errors(caplog)[-1]


def test_a_report_that_cannot_be_written_exits_one(written, dataset, tmp_path, caplog):
    outdir, _scores = written
    assert audit_main(outdir, reference=dataset.reference, report=str(tmp_path / "no" / "r.md")) == 1
    assert "r.md" in errors(caplog)[-1]
    assert not (tmp_path / "no").exists()
