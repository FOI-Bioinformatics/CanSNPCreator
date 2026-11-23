"""Merging schemes built in tiers: a scheme over the whole group whose leaf the
finer scheme opens up, joined at that one shared node, the seam.

The tests split a written scheme at a node and merge the two halves back. A merge
that is right gives back the scheme it started from; each refusal is a way two
halves can fail to fit.
"""

import logging
import os

import pytest
from conftest import METADATA, cli_run
from test_export import needs_validator, run_validator, scheme_checkout

from cansnpscorer.scheme import (
    CANSNP_COLUMNS,
    CANSNP_SUFFIX,
    GENOMES_SUFFIX,
    REFERENCES_SUFFIX,
    TREE_SUFFIX,
    read_scheme,
)

ORGANISM = METADATA["organism"]


def lines(path):
    with open(path) as handle:
        return handle.read().splitlines()


def write(path, rows):
    with open(path, "w", newline="\n") as handle:
        handle.writelines(row + "\n" for row in rows)


def split_at(outdir, seam, tmp_path):
    """Split a written scheme at ``seam`` into a coarse scheme, where the seam is
    a leaf, and a fine scheme rooted at it. Both keep the seam's canSNPs."""
    tree = [line.split("\t") for line in lines(os.path.join(outdir, ORGANISM + TREE_SUFFIX))]
    below = {node for path in tree if seam in path for node in path[path.index(seam) + 1:]}
    coarse_paths = [path for path in tree if not below & set(path)]
    fine_paths = [path[path.index(seam):] for path in tree if seam in path]
    coarse_nodes = {node for path in coarse_paths for node in path}
    fine_nodes = {node for path in fine_paths for node in path}

    cansnps = lines(os.path.join(outdir, ORGANISM + CANSNP_SUFFIX))
    header, rows = cansnps[0], cansnps[1:]
    references = lines(os.path.join(outdir, ORGANISM + REFERENCES_SUFFIX))

    parts = {}
    for name, paths, nodes in (("coarse", coarse_paths, coarse_nodes),
                               ("fine", fine_paths, fine_nodes)):
        directory = tmp_path / name
        directory.mkdir()
        write(directory / (ORGANISM + TREE_SUFFIX), ["\t".join(path) for path in paths])
        write(directory / (ORGANISM + CANSNP_SUFFIX),
              [header, *[row for row in rows if row.split("\t")[0] in nodes]])
        write(directory / (ORGANISM + REFERENCES_SUFFIX), references)
        parts[name] = str(directory)
    return parts["coarse"], parts["fine"]


def merge(tmp_path, *schemes, outdir=None):
    outdir = outdir or str(tmp_path / "merged")
    arguments = ["merge"]
    for scheme in schemes:
        arguments += ["--scheme", scheme]
    return cli_run(*arguments, "--outdir", outdir, "--organism", ORGANISM), outdir


def cansnp_rows(directory):
    """The canSNPs of a scheme, without the file's own line numbers."""
    return sorted((s.node, s.rank, s.sequence, s.position, s.derived, s.ancestral, s.typed)
                  for s in read_scheme(directory).snp_set)


def edit(path, replace, with_text):
    with open(path) as handle:
        text = handle.read()
    assert replace in text, replace
    with open(path, "w") as handle:
        handle.write(text.replace(replace, with_text, 1))


def cansnp_line(directory, node):
    """The first canSNP line of a node, split into its columns."""
    for line in lines(os.path.join(directory, ORGANISM + CANSNP_SUFFIX))[1:]:
        fields = line.split("\t")
        if fields[0] == node:
            return dict(zip(CANSNP_COLUMNS, fields, strict=True))
    raise AssertionError(node)


# -- the round trip -------------------------------------------------------------


def test_a_scheme_split_at_a_node_merges_back_to_itself(written, tmp_path):
    outdir, _scores = written
    coarse, fine = split_at(outdir, "A", tmp_path)
    status, merged = merge(tmp_path, coarse, fine)
    assert status == 0
    assert cansnp_rows(merged) == cansnp_rows(outdir)
    original = read_scheme(outdir)
    result = read_scheme(merged)
    assert result.problems == []
    assert set(result.nodes) == set(original.nodes)
    assert {node: result.hierarchy.parents[node] for node in result.nodes} == \
        {node: original.hierarchy.parents[node] for node in original.nodes}


def test_the_references_are_written_once(written, tmp_path):
    outdir, _scores = written
    coarse, fine = split_at(outdir, "A", tmp_path)
    _status, merged = merge(tmp_path, coarse, fine)
    assert lines(os.path.join(merged, ORGANISM + REFERENCES_SUFFIX)) == \
        lines(os.path.join(outdir, ORGANISM + REFERENCES_SUFFIX))


def test_the_seam_keeps_the_canSNPs_of_the_first_scheme(written, tmp_path, caplog):
    # Only the coarse scheme has a genome outside the group to settle the seam's
    # polarity, so its canSNPs win; the fine scheme's are dropped, and that is said.
    caplog.set_level(logging.INFO)
    outdir, _scores = written
    coarse, fine = split_at(outdir, "A", tmp_path)
    seam = cansnp_line(fine, "A")
    position = int(seam["sequence_position"])
    edit(os.path.join(fine, ORGANISM + CANSNP_SUFFIX),
         f"\t{seam['sequence']}\t{position}\t", f"\t{seam['sequence']}\t{position + 3}\t")
    status, merged = merge(tmp_path, coarse, fine)
    assert status == 0
    assert cansnp_line(merged, "A")["sequence_position"] == str(position)
    assert "seam" in caplog.text and "dropped" in caplog.text


@needs_validator
def test_the_merged_scheme_passes_the_cansnper2_data_validator(written, tmp_path):
    outdir, _scores = written
    coarse, fine = split_at(outdir, "A", tmp_path)
    write(os.path.join(coarse, ORGANISM + GENOMES_SUFFIX),
          ["genome\tnode\trole", "GCF_000000001.1\tA\trepresentative",
           "GCF_000000009.1\tB.1\trepresentative"])
    write(os.path.join(fine, ORGANISM + GENOMES_SUFFIX),
          ["genome\tnode\trole", "GCF_000000001.1\tA.1\trepresentative"])
    status, merged = merge(tmp_path, coarse, fine)
    assert status == 0
    completed = run_validator(scheme_checkout(merged, tmp_path))
    assert completed.returncode == 0, completed.stdout + completed.stderr


# -- what is refused -------------------------------------------------------------


def test_a_site_used_by_two_nodes_is_refused(written, tmp_path, caplog):
    # The tiers were extracted independently from one reference, so two of them
    # can pick one site for different nodes; CanSNPer2 refuses a repeated site.
    outdir, _scores = written
    coarse, fine = split_at(outdir, "A", tmp_path)
    taken = cansnp_line(coarse, "B")
    mine = cansnp_line(fine, "A.1")
    edit(os.path.join(fine, ORGANISM + CANSNP_SUFFIX),
         f"\t{mine['sequence']}\t{mine['sequence_position']}\t",
         f"\t{taken['sequence']}\t{taken['sequence_position']}\t")
    status, merged = merge(tmp_path, coarse, fine)
    assert status == 2
    assert not os.path.exists(merged), "nothing is written"
    assert f"{taken['sequence']}:{taken['sequence_position']}" in caplog.text


def test_schemes_that_share_more_than_the_seam_are_refused(written, tmp_path, caplog):
    outdir, _scores = written
    coarse, fine = split_at(outdir, "A", tmp_path)
    with open(os.path.join(fine, ORGANISM + TREE_SUFFIX), "a") as handle:
        handle.write("A\tB.1\n")
    status, _merged = merge(tmp_path, coarse, fine)
    assert status == 2
    assert "B.1" in caplog.text


def test_a_scheme_that_shares_no_node_is_refused(written, tmp_path, caplog):
    # Its tree would hang from nothing, and the merged tree would have two roots.
    outdir, _scores = written
    coarse, fine = split_at(outdir, "A", tmp_path)
    for suffix in (TREE_SUFFIX, CANSNP_SUFFIX):
        path = os.path.join(fine, ORGANISM + suffix)
        write(path, [line.replace("A", "Z", 1) if line.startswith("A") else line
                     for line in lines(path)])
    status, _merged = merge(tmp_path, coarse, fine)
    assert status == 2
    assert "shares no node" in caplog.text


def test_references_that_disagree_are_refused(written, tmp_path, caplog):
    outdir, _scores = written
    coarse, fine = split_at(outdir, "A", tmp_path)
    edit(os.path.join(fine, ORGANISM + REFERENCES_SUFFIX), "ASM712v1", "ASM999v1")
    status, _merged = merge(tmp_path, coarse, fine)
    assert status == 2
    assert "reference" in caplog.text


def test_a_scheme_with_problems_of_its_own_is_refused(written, tmp_path, caplog):
    outdir, _scores = written
    coarse, fine = split_at(outdir, "A", tmp_path)
    row = cansnp_line(fine, "A.1")
    edit(os.path.join(fine, ORGANISM + CANSNP_SUFFIX), f"A.1\t{row['rank']}\t", "A.1\t7\t")
    status, _merged = merge(tmp_path, coarse, fine)
    assert status == 2
    assert "A.1" in caplog.text


# -- the genomes file -----------------------------------------------------------


def test_genomes_at_the_seam_give_way_to_the_finer_scheme(written, tmp_path):
    # In the coarse scheme the seam is a leaf, so its genomes file lists the fine
    # scheme's genomes there. After the merge the seam is an inner node, and the
    # genomes file lists terminal nodes only.
    outdir, _scores = written
    coarse, fine = split_at(outdir, "A", tmp_path)
    write(os.path.join(coarse, ORGANISM + GENOMES_SUFFIX),
          ["genome\tnode\trole", "GCF_000000001.1\tA\trepresentative",
           "GCF_000000009.1\tB.1\trepresentative", "GCF_000000099.1\t\toutgroup"])
    write(os.path.join(fine, ORGANISM + GENOMES_SUFFIX),
          ["genome\tnode\trole", "GCF_000000001.1\tA.1\trepresentative",
           "GCF_000000009.1\t\toutgroup"])
    status, merged = merge(tmp_path, coarse, fine)
    assert status == 0
    genomes = sorted(read_scheme(merged).genomes)
    assert genomes == [
        ("GCF_000000001.1", "A.1", "representative"),
        ("GCF_000000009.1", "B.1", "representative"),
        ("GCF_000000099.1", "", "outgroup"),
    ], "the fine scheme's outgroup is placed by the coarse one, and stays placed"


def test_a_genome_at_two_terminal_nodes_is_refused(written, tmp_path, caplog):
    outdir, _scores = written
    coarse, fine = split_at(outdir, "A", tmp_path)
    write(os.path.join(coarse, ORGANISM + GENOMES_SUFFIX),
          ["genome\tnode\trole", "GCF_000000001.1\tB.1\trepresentative"])
    write(os.path.join(fine, ORGANISM + GENOMES_SUFFIX),
          ["genome\tnode\trole", "GCF_000000001.1\tA.1\trepresentative"])
    status, _merged = merge(tmp_path, coarse, fine)
    assert status == 2
    assert "GCF_000000001.1" in caplog.text


@pytest.mark.parametrize("missing", ["nowhere"])
def test_a_scheme_folder_that_does_not_exist_exits_one(tmp_path, missing):
    status, _merged = merge(tmp_path, str(tmp_path / missing), str(tmp_path / missing))
    assert status == 1
