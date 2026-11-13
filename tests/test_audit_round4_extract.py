"""Defects of extract and evaluate found in the fourth audit round, each pinned
by the input that showed it."""

import logging
import random

import pytest
from conftest import cli_run, vcf_lines, vcf_records

from cansnpscorer import cansnp_extractor, vcf
from cansnpscorer.helpers.data_utils import gc_content, reference_anchors
from cansnpscorer.helpers.file_utils import read_hierarchy

TWO_SIDES = "child\tparent\nA\troot\nB\troot\n"


def random_sequence(length, seed=1):
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(length))


def changed(sequence, *positions):
    other = {"A": "T", "T": "A", "C": "G", "G": "C"}
    bases = list(sequence)
    for position in positions:
        bases[position] = other[bases[position]]
    return "".join(bases)


def extract(tmp_path, reference, alignment, genome_map, hierarchy=TWO_SIDES):
    """Run extract on these records; the output path and its records."""
    for name, records in (("ref.fa", reference), ("aln.fa", alignment)):
        (tmp_path / name).write_text("".join(f">{n}\n{s}\n" for n, s in records))
    (tmp_path / "h.tsv").write_text(hierarchy)
    (tmp_path / "g.tsv").write_text(genome_map)
    output = tmp_path / "out.vcf"
    status = cansnp_extractor.main(str(tmp_path / "aln.fa"), str(tmp_path / "ref.fa"),
                                   str(tmp_path / "h.tsv"), str(tmp_path / "g.tsv"),
                                   str(output), 1.0)
    return status, output


# -- extract -----------------------------------------------------------------


@pytest.mark.parametrize("layout", ["insertion", "deletion"])
def test_a_snp_beside_a_variable_indel_is_not_taken_as_isolated(tmp_path, layout):
    sequence = random_sequence(60)
    derived = changed(sequence, 10)
    genome_map = "a1\tA\na2\tA\nb1\tB\nb2\tB\n"
    if layout == "insertion":
        # Three columns after reference position 12 that only the other genomes
        # fill, and with different bases: a gapped layout.
        def insert(row, text):
            return row[:12] + text + row[12:]
        alignment = [("REF", insert(sequence, "---")), ("og", insert(sequence, "---")),
                     ("a1", insert(derived, "ACG")), ("a2", insert(derived, "ATG")),
                     ("b1", insert(sequence, "ACG")), ("b2", insert(sequence, "AGG"))]
    else:
        # Three bases deleted in one genome of B, right beside the canSNP.
        alignment = [("REF", sequence), ("og", sequence), ("a1", derived), ("a2", derived),
                     ("b1", sequence), ("b2", sequence[:11] + "---" + sequence[14:])]
    status, output = extract(tmp_path, [("c", sequence)], alignment, genome_map)
    assert status == 0
    (record,) = [r for r in vcf_records(output) if r.node == "A"]
    assert record.pos == 11 and record.info["NEAREST"] == "1"


def test_a_gapped_column_counts_at_the_reference_bases_on_either_side():
    # Columns 3 to 6 are an insertion: reference bases 0 1 2 | gap x4 | 3 4 5.
    positions = [0, 1, 2, None, None, None, None, 3, 4, 5]
    anchors = reference_anchors([4, 5, 9], positions.__getitem__, len(positions))
    assert anchors == [2, 3, 5]
    # An insertion at either end of the alignment has one side only.
    positions = [None, 0, 1, None]
    assert reference_anchors([0, 3], positions.__getitem__, len(positions)) == [0, 1]


def test_an_illegal_contig_name_is_refused_before_any_work(tmp_path, caplog):
    sequence = random_sequence(60)
    alignment = [("REF", sequence), ("a", changed(sequence, 5)), ("b", sequence)]
    with caplog.at_level(logging.INFO):
        status, output = extract(tmp_path, [("contig(1)", sequence)], alignment,
                                 "a\tA\nb\tB\n")
    assert status == 1 and not output.exists()
    assert any("does not allow" in m and "contig(1)" in m for m in caplog.messages)
    assert not any("SNP extraction" in m for m in caplog.messages), "refused before the work"


def test_spaces_around_the_hierarchy_header_cells_are_ignored(tmp_path):
    path = tmp_path / "h.tsv"
    path.write_text(" child\tparent \nA\troot\nB\troot\nA1\tA\n")
    hierarchy = read_hierarchy(str(path))
    assert hierarchy.root == "root" and hierarchy.parents["A1"] == "A"


def test_root_polarity_does_not_follow_the_order_of_rows(tmp_path):
    # The reference sits in A and no genome lies outside A and B, so the column
    # cannot tell which side is derived: the reference's base is taken as
    # ancestral, and the SNP marks B whichever side the hierarchy lists first.
    sequence = random_sequence(60)
    derived = changed(sequence, 10)
    alignment = [("REF", sequence), ("a1", sequence), ("b1", derived), ("b2", derived)]
    found = []
    for number, hierarchy in enumerate((TWO_SIDES, "child\tparent\nB\troot\nA\troot\n")):
        directory = tmp_path / str(number)
        directory.mkdir()
        status, output = extract(directory, [("c", sequence)], alignment,
                                 "REF\tA\na1\tA\nb1\tB\nb2\tB\n", hierarchy)
        assert status == 0
        found.append([(r.pos, r.node, r.ancestral) for r in vcf_records(output) if r.node])
    assert found[0] == found[1] == [(11, "B", sequence[10])]


# -- evaluate ----------------------------------------------------------------


def evaluate(dataset, tmp_path, *extra):
    snps = str(tmp_path / "snps.vcf")
    assert cli_run("extract", "--genome-map", dataset.genome_map, "--hierarchy",
                   dataset.hierarchy, "--alignment", dataset.alignment, "--reference",
                   dataset.reference, "-o", snps) == 0
    scores = str(tmp_path / "scores.vcf")
    assert cli_run("evaluate", "--snps", snps, "--reference", dataset.reference,
                   "--genes", dataset.gff, "--repeats", dataset.repeats, "-o", scores,
                   *extra) == 0
    return scores


def test_a_flank_beyond_the_nearest_cap_is_warned_of(dataset, tmp_path, caplog):
    with caplog.at_level(logging.WARNING):
        evaluate(dataset, tmp_path, "--flank", str(vcf.NEAREST_CAP))
    assert not any("--flank" in m for m in caplog.messages)
    caplog.clear()
    with caplog.at_level(logging.WARNING):
        evaluate(dataset, tmp_path, "--flank", str(2 * vcf.NEAREST_CAP))
    assert any("--flank" in m and "not known" in m for m in caplog.messages)


def test_gc_is_taken_over_a_c_g_and_t_only():
    window = "N" * 45 + "GCGCGCATAT" + "N" * 45
    assert gc_content(50, window, 0.6) == pytest.approx(gc_content(5, "GCGCGCATAT", 0.6))
    assert gc_content(50, window, 0.6) == pytest.approx(1.0), "6 of 10 called bases are GC"
    assert gc_content(50, "N" * 100, 0.5) == 0.0, "nothing known of the window"


def test_the_nottop_description_names_how_the_set_is_chosen(dataset, tmp_path):
    scores = evaluate(dataset, tmp_path, "--top-n", "1", "--min-spacing", "250")
    header, _records = vcf_lines(scores)
    (line,) = [h for h in header if h.startswith("##FILTER=<ID=NotTop")]
    assert "reliability" in line and "250 bases" in line and "different genes" in line
    unspaced = next(f for f in vcf.score_filters(0, 1, 0) if f.id == "NotTop")
    assert "bases of one already chosen" not in unspaced.description
