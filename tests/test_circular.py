"""Features across the origin of a circular contig, as NCBI's GFF3 gives them."""

import logging

from conftest import cli_run

from cansnpscorer.cansnp_evaluator import wrap_origin

CODON = 3


def test_a_gene_across_the_origin_is_split_in_two():
    kept, split = wrap_origin([(250, 319, "+")], 300)
    assert split == 1
    assert kept == [(0, 19, "+"), (250, 299, "+")]


def test_a_feature_within_the_contig_is_left_alone():
    kept, split = wrap_origin([(10, 99, "+")], 300)
    assert (kept, split) == ([(10, 99, "+")], 0)


def test_a_feature_longer_than_the_contig_is_not_split():
    # Its two parts would overlap: it is malformed, not circular.
    kept, split = wrap_origin([(10, 4999, "+")], 300)
    assert (kept, split) == ([(10, 4999, "+")], 0)


def codon_starts(parts, strand):
    """The positions where codons begin, along the direction of reading.

    Across the origin, a + strand CDS is read from the part at the contig's
    end and then from its start; a - strand CDS the other way round."""
    order = sorted(parts, key=lambda part: part[0], reverse=strand == "+")
    starts = []
    for start, end, _strand, phase in order:
        bases = range(start, end + 1) if strand == "+" else range(end, start - 1, -1)
        starts += list(bases)[phase::CODON]
    return starts


def test_a_cds_across_the_origin_keeps_its_codons():
    length = 10
    for strand in "+-":
        for first_base in (6, 7, 8):
            for phase in (0, 1, 2):
                end = length + 4
                parts, _split = wrap_origin([(first_base, end, strand, phase)], length,
                                            cds=True)
                # The bases of the whole CDS in reading order, around the origin.
                bases = [position % length for position in range(first_base, end + 1)]
                if strand == "-":
                    bases.reverse()
                expected = bases[phase::CODON]
                assert codon_starts(parts, strand) == expected, (strand, first_base, phase)


def test_evaluate_reads_a_gff_with_a_gene_across_the_origin(written, dataset, tmp_path,
                                                             caplog):
    gff = tmp_path / "circular.gff"
    with open(dataset.gff) as handle:
        lines = handle.read()
    # c2 is 300 bases long: a gene and its CDS from 251 over the origin to 20.
    lines += ("c2\tRefSeq\tgene\t251\t320\t.\t+\t.\tID=gene-wrap\n"
              "c2\tRefSeq\tCDS\t251\t320\t.\t+\t0\tID=cds-wrap;Parent=gene-wrap\n")
    gff.write_text(lines)
    with caplog.at_level(logging.WARNING):
        assert cli_run("evaluate", "--snps", str(tmp_path / "snps.vcf"), "--reference",
                       dataset.reference, "--genes", str(gff), "--repeats", dataset.repeats,
                       "-o", str(tmp_path / "wrapped.vcf")) == 0
    assert "across the origin of a circular sequence" in caplog.text


def test_evaluate_still_refuses_a_feature_longer_than_its_contig(written, dataset, tmp_path):
    gff = tmp_path / "broken.gff"
    with open(dataset.gff) as handle:
        original = handle.read()
    gff.write_text(original + "c2\tRefSeq\tgene\t11\t5000\t.\t+\t.\tID=gene-long\n")
    assert cli_run("evaluate", "--snps", str(tmp_path / "snps.vcf"), "--reference",
                   dataset.reference, "--genes", str(gff), "--repeats", dataset.repeats,
                   "-o", str(tmp_path / "broken.vcf")) == 2
