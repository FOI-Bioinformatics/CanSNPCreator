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

SNP_OLD = "testus_organismus_snp.txt"


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
    order = [4, 5, 0, 1, 2, 3, 6, 7, 8, 9]
    with open(path, "w") as handle:
        handle.writelines("\t".join(row[index] for index in order) + "\n" for row in rows)
    scheme = read_scheme(outdir)
    assert not scheme.problems
    assert {snp.node: (snp.sequence, snp.position, snp.derived) for snp in scheme.snps}["A"] \
        == ("c1", 201, "C")


def old_snp_file(outdir, header):
    """Replace the canSNP file by a _snp.txt, as an earlier export wrote it."""
    os.remove(os.path.join(outdir, SNPS))
    path = os.path.join(outdir, SNP_OLD)
    with open(path, "w") as handle:
        handle.write(header + "\n")
    return path


def test_old_header_names_the_conversion_command(written):
    outdir, _scores = written
    path = old_snp_file(
        outdir, "snp_id\tstrain\treference\tgenome\tposition\tderived_base\tancestral_base")
    with open(path, "a") as handle:
        handle.write("A\ttestus\tx\tTEST1\t201\tC\tG\n")
    with pytest.raises(SchemeError, match=r"CanSNPer2-database --convert") as raised:
        read_scheme(outdir)
    assert SNP_OLD in str(raised.value)
    assert "2.3.0" in str(raised.value)


def test_old_snp_set_header_names_the_conversion_command(written):
    outdir, _scores = written
    with open(os.path.join(outdir, SNPS), "w") as handle:
        handle.write("node\trank\tgenome\tposition\tderived_base\tancestral_base\tscore\n"
                     "A\t1\tTEST1\t201\tC\tG\t500\n")
    with pytest.raises(SchemeError, match=r"CanSNPer2-database --convert") as raised:
        read_scheme(outdir)
    assert SNPS in str(raised.value)


def test_a_snp_file_of_an_earlier_export_is_read_as_rank_one_and_typed(written):
    outdir, _scores = written
    path = old_snp_file(
        outdir, "snp_id\tstrain\treference\tgenome\tsequence\tsequence_position\t"
                "derived_base\tancestral_base")
    with open(path, "a") as handle:
        handle.write("A\ttestus\tSjodin2026\tTEST1\tc1\t201\tC\tG\n")
    scheme = read_scheme(outdir)
    assert not [problem for problem in scheme.problems if "ranks" in problem or "typed" in problem]
    [snp] = scheme.snps
    assert (snp.node, snp.rank, snp.typed, snp.score, snp.position) == ("A", 1, True, "", 201)
    assert scheme.snp_set == scheme.snps


def test_a_missing_column_is_a_problem(written):
    outdir, _scores = written
    edit(outdir, SNPS, "\tsequence\t", "\tcontig\t")
    scheme = read_scheme(outdir)
    assert not scheme.snps
    assert any("sequence" in problem for problem in scheme.problems)


def rewrite_typed(outdir, typed_by_row):
    """Set the typed column of the canSNP file, row by row in file order."""
    path = os.path.join(outdir, SNPS)
    with open(path) as handle:
        lines = handle.read().splitlines()
    rows = [line.split("\t") for line in lines[1:]]
    for row, value in zip(rows, typed_by_row, strict=True):
        row[-1] = value
    with open(path, "w") as handle:
        handle.write("\n".join([lines[0], *("\t".join(row) for row in rows)]) + "\n")


def test_a_node_with_no_typed_cansnp_is_a_problem(written):
    outdir, _scores = written
    with open(os.path.join(outdir, SNPS)) as handle:
        count = len(handle.read().splitlines()) - 1
    rewrite_typed(outdir, ["no"] + ["yes"] * (count - 1))
    scheme = read_scheme(outdir)
    assert any("marks no canSNP of node 'root' typed" in problem for problem in scheme.problems)
    assert "root" not in scheme.snps_by_node
    assert any(snp.node == "root" for snp in scheme.snp_set)


def test_a_node_with_two_typed_cansnps_is_a_problem(written):
    outdir, _scores = written
    with open(os.path.join(outdir, SNPS), "a") as handle:
        handle.write("A\t2\tSjodin2026\tTEST1\tc1\t211\tC\tG\t100\tyes\n")
    scheme = read_scheme(outdir)
    assert any("marks 2 canSNPs of node 'A' typed" in problem for problem in scheme.problems)
    assert not any("'root'" in problem for problem in scheme.problems)


@pytest.mark.parametrize("value", ["maybe", "1", ""])
def test_a_typed_value_that_is_not_yes_or_no_is_a_problem(written, value):
    outdir, _scores = written
    edit(outdir, SNPS, "A.1\t1\tSjodin2026\tTEST1\tc1\t401\tG\tC\t350\tyes",
         f"A.1\t1\tSjodin2026\tTEST1\tc1\t401\tG\tC\t350\t{value}")
    scheme = read_scheme(outdir)
    assert any("line 4" in problem and "where yes or no is expected" in problem
               for problem in scheme.problems)
    assert "A.1" not in {snp.node for snp in scheme.snp_set}


def test_the_typed_cansnp_need_not_be_rank_one(written):
    outdir, _scores = written
    edit(outdir, SNPS, "A\t1\tSjodin2026\tTEST1\tc1\t201\tC\tG\t350\tyes",
         "A\t1\tSjodin2026\tTEST1\tc1\t201\tC\tG\t350\tno\n"
         "A\t2\tSjodin2026\tTEST1\tc1\t211\tC\tG\t100\tyes")
    scheme = read_scheme(outdir)
    assert not scheme.problems
    assert [(snp.rank, snp.position) for snp in scheme.snps if snp.node == "A"] == [(2, 211)]
    assert [snp.rank for snp in scheme.snp_set if snp.node == "A"] == [1, 2]
    [alternative] = [snp for snp in scheme.snp_set if snp.node == "A" and not snp.typed]
    assert alternative.tag("x").startswith("the canSNP of rank 1: ")
