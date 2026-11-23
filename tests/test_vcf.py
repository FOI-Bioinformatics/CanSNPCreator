"""VCF 4.3 as this tool writes and reads it."""

from pathlib import Path

import pytest

from cansnpscorer import vcf
from cansnpscorer.model import Reference

REFERENCE = Reference.from_records([("chr1", "ACGTACGTAC"), ("pX", "GGGCCC")])


def written(tmp_path, records, infos=(vcf.NODE, vcf.ANCESTRAL)):
    path = str(tmp_path / "snps.vcf")
    vcf.write(path, vcf.header_lines(REFERENCE, "ref.fasta", "test", list(infos)), records)
    return path


def test_ref_is_the_reference_base_and_aa_gives_the_polarity(tmp_path):
    # The reference carries the derived base at chr1:2, so REF is the derived
    # allele and AA names the other.
    record = vcf.snp_record(REFERENCE, "chr1", 2, ("C", "T"), node="A.1", ancestral="T")
    assert (record.ref, record.alt, record.ancestral, record.derived) == ("C", "T", "T", "C")
    path = written(tmp_path, [record])
    read = vcf.read(path).records[0]
    assert (read.chrom, read.pos, read.node, read.derived, read.ancestral) == ("chr1", 2, "A.1", "C", "T")


def test_a_snp_whose_alleles_miss_the_reference_base_is_refused():
    with pytest.raises(vcf.VcfError, match="neither allele"):
        vcf.snp_record(REFERENCE, "chr1", 1, ("C", "G"))


def test_reserved_characters_in_a_node_name_survive(tmp_path):
    node = "A;B=C,D%E:F"
    record = vcf.snp_record(REFERENCE, "chr1", 1, ("A", "G"), node=node, ancestral="A")
    path = written(tmp_path, [record])
    assert "NODE=A%3BB%3DC%2CD%25E%3AF" in Path(path).read_text()
    assert vcf.read(path).records[0].node == node


def test_the_header_declares_contigs_and_fields(tmp_path):
    path = written(tmp_path, [])
    text = Path(path).read_text()
    assert text.startswith("##fileformat=VCFv4.3\n")
    assert "##contig=<ID=chr1,length=10>\n##contig=<ID=pX,length=6>" in text
    assert '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">' in text
    parsed = vcf.read(path)
    assert parsed.contigs == {"chr1": 10, "pX": 6}
    assert set(parsed.infos) == {"NODE", "AA"}


@pytest.mark.parametrize(
    "record,message",
    [
        ("chr1\t1\t.\tA\tG,T\t.\tPASS\t.", "not a biallelic SNP"),
        ("chr1\t1\t.\tAC\tG\t.\tPASS\t.", "not a biallelic SNP"),
        ("chr1\t1\t.\tA\tG\t.\tPASS\tDP=3", "does not declare"),
        ("chr1\tone\t.\tA\tG\t.\tPASS\t.", "position"),
    ],
)
def test_what_this_tool_did_not_write_is_refused(tmp_path, record, message):
    path = written(tmp_path, [])
    with open(path, "a") as handle:
        handle.write(record + "\n")
    with pytest.raises(vcf.VcfError, match=message):
        vcf.read(path)


def test_records_are_checked_against_the_reference(tmp_path):
    path = written(tmp_path, [])
    with open(path, "a") as handle:
        handle.write("chr1\t1\t.\tG\tT\t.\tPASS\t.\n")
    with pytest.raises(vcf.VcfError, match="has REF G where the reference has A"):
        vcf.read(path).check_against(REFERENCE, path)

    other = Reference.from_records([("chr1", "ACGTACGTACG"), ("pX", "GGGCCC")])
    good = written(tmp_path, [])
    with pytest.raises(vcf.VcfError, match="another reference"):
        vcf.read(good).check_against(other, good)


def test_a_file_that_is_not_vcf_is_refused(tmp_path):
    path = tmp_path / "x.vcf"
    path.write_text("contig\tposition\n")
    with pytest.raises(vcf.VcfError, match="not VCF"):
        vcf.read(str(path))
