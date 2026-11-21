"""Reading a CanSNPer2-data scheme in the format of CanSNPer2 2.3.0.

From CanSNPer2 2.3.0 a canSNP is given as a reference sequence and a 1-based
position on that sequence. Files written before then give one position along
the reference sequences joined in file order, and are refused with the command
that converts them.
"""

import os

import pytest
from conftest import PLANTED, SNPS, edit

from cansnpscorer.scheme import SchemeError, read_scheme

SNP_SET = "testus_organismus_canSNPs.txt"


def test_snps_are_read_by_sequence_and_position(written, dataset):
    outdir, _scores = written
    scheme = read_scheme(outdir)
    assert not scheme.problems
    found = {snp.node: (snp.sequence, snp.position) for snp in scheme.snps}
    assert found == {node: dataset.contig_of(column)
                     for node, (column, _ancestral, _derived) in PLANTED.items()}
    assert found["B.2"] == ("c2", 1), "the first base of the second sequence"


def test_columns_are_read_by_name(written):
    # A file whose columns come in another order is read by the header.
    outdir, _scores = written
    path = os.path.join(outdir, SNPS)
    with open(path) as handle:
        rows = [line.split("\t") for line in handle.read().splitlines()]
    order = [4, 5, 0, 1, 2, 3, 6, 7]
    with open(path, "w") as handle:
        handle.writelines("\t".join(row[index] for index in order) + "\n" for row in rows)
    scheme = read_scheme(outdir)
    assert not scheme.problems
    assert {snp.node: (snp.sequence, snp.position, snp.derived) for snp in scheme.snps}["A"] \
        == ("c1", 201, "C")


def test_old_header_names_the_conversion_command(written):
    outdir, _scores = written
    path = os.path.join(outdir, SNPS)
    with open(path) as handle:
        lines = handle.read().splitlines()
    old = ["snp_id\tstrain\treference\tgenome\tposition\tderived_base\tancestral_base"]
    for line in lines[1:]:
        fields = line.split("\t")
        old.append("\t".join(fields[:4] + fields[5:]))
    with open(path, "w") as handle:
        handle.write("\n".join(old) + "\n")

    with pytest.raises(SchemeError, match=r"CanSNPer2-database --convert") as raised:
        read_scheme(outdir)
    assert SNPS in str(raised.value)
    assert "2.3.0" in str(raised.value)


def test_old_snp_set_header_names_the_conversion_command(written):
    outdir, _scores = written
    with open(os.path.join(outdir, SNP_SET), "w") as handle:
        handle.write("node\trank\tgenome\tposition\tderived_base\tancestral_base\tscore\n"
                     "A\t1\tTEST1\t201\tC\tG\t500\n")
    with pytest.raises(SchemeError, match=r"CanSNPer2-database --convert") as raised:
        read_scheme(outdir)
    assert SNP_SET in str(raised.value)


def test_a_missing_column_is_a_problem(written):
    outdir, _scores = written
    edit(outdir, SNPS, "\tsequence\t", "\tcontig\t")
    scheme = read_scheme(outdir)
    assert not scheme.snps
    assert any("sequence" in problem for problem in scheme.problems)


def test_the_snp_set_rank_one_must_match_on_the_sequence(written):
    outdir, _scores = written
    scheme = read_scheme(outdir)
    a = next(snp for snp in scheme.snps if snp.node == "A")
    with open(os.path.join(outdir, SNP_SET), "w") as handle:
        handle.write("node\trank\tgenome\tsequence\tsequence_position\tderived_base\t"
                     "ancestral_base\tscore\n"
                     f"A\t1\t{a.genome}\tc2\t{a.position}\t{a.derived}\t{a.ancestral}\t500\n")
    scheme = read_scheme(outdir)
    assert any("rank 1 that is not the SNP" in problem for problem in scheme.problems)
