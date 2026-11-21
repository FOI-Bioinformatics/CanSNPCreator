"""The audit against the genomes themselves.

The scheme is the one the export subcommand writes from the planted dataset, so
every SNP in it is canonical in that dataset's alignment. Each test changes the
alignment, the group file or the scheme in one way and asserts what the audit
says about it.
"""

import os

from conftest import (
    CLUSTERS,
    INCONGRUENT,
    PLANTED,
    SNPS,
    cli_run,
    edit,
    vcf_records,
    write_fasta,
)

from cansnpscorer.cansnp_audit import PROBLEM, audit_scheme
from cansnpscorer.cansnp_audit import main as audit_main
from cansnpscorer.helpers.file_utils import (
    read_alignment,
    read_genome_map,
    reference_parser,
)
from cansnpscorer.scheme import read_scheme


def audit(outdir, dataset, alignment_path=None, group_path=None, **arguments):
    alignment = read_alignment(alignment_path or dataset.alignment)
    group_path = group_path or dataset.genome_map
    clusters = read_genome_map(group_path)
    return audit_scheme(read_scheme(outdir), reference_parser(dataset.reference),
                        alignment=alignment, clusters=clusters, **arguments)


def of(findings, category, severity=None):
    return [finding for finding in findings if finding.category == category
            and (severity is None or finding.severity == severity)]


def changed_alignment(dataset, tmp_path, changes):
    """Write the alignment with some bases changed, as {genome: {column: base}}."""
    records = read_alignment(dataset.alignment)
    for name, columns in changes.items():
        sequence = list(records[name])
        for column, base in columns.items():
            sequence[column] = base
        records[name] = "".join(sequence)
    path = str(tmp_path / "changed.fasta")
    write_fasta(path, list(records.items()))
    return path


def write_groups(tmp_path, placements):
    path = tmp_path / "genome_map.tsv"
    path.write_text("".join(
        f"{genome}\t{node}\n" for node, genomes in placements.items() for genome in genomes))
    return str(path)


# -- a sound scheme ----------------------------------------------------------


def test_every_planted_snp_is_canonical_in_the_genomes(written, dataset):
    outdir, _scores = written
    findings = audit(outdir, dataset)
    for category in ("not-canonical", "homoplasy", "polarity", "unsupported-node"):
        assert of(findings, category) == [], category
    assert [finding.severity for finding in findings if finding.severity == PROBLEM] == []


def test_the_reference_is_reported_as_placed_nowhere(written, dataset):
    outdir, _scores = written
    ungrouped = of(audit(outdir, dataset), "ungrouped")
    assert len(ungrouped) == 1 and "ref" in ungrouped[0].detail


def test_a_node_marked_by_one_site_only_is_reported(written, dataset):
    # The planted data carries exactly one canSNP per node, so every branch
    # rests on a single base.
    outdir, _scores = written
    single = {finding.node for finding in of(audit(outdir, dataset), "single-site")}
    assert single == set(PLANTED) - {"root"}


def test_the_incongruent_site_is_counted_as_conflicting(written, dataset):
    outdir, _scores = written
    conflicting = of(audit(outdir, dataset), "conflicting-sites")
    assert len(conflicting) == 1
    assert conflicting[0].detail.startswith("1 of ")


# -- a SNP that is not canonical ---------------------------------------------


def test_an_isolate_without_the_derived_base_is_reported(written, dataset, tmp_path):
    column, ancestral, _derived = PLANTED["A.1"]
    path = changed_alignment(dataset, tmp_path, {"s1b": {column: ancestral}})
    findings = audit(written[0], dataset, alignment_path=path)

    failing = of(findings, "not-canonical", PROBLEM)
    assert [finding.node for finding in failing] == ["A.1"]
    assert "1 of 2 isolates" in failing[0].detail and "s1b" in failing[0].detail
    # Nothing else in the alignment separates s1 and s1b from the rest, so the
    # node is no longer a clade in these genomes.
    assert [finding.node for finding in of(findings, "unsupported-node")] == ["A.1"]


def test_a_genome_outside_the_node_with_the_derived_base_is_reported(written, dataset,
                                                                      tmp_path):
    column, _ancestral, derived = PLANTED["A.1"]
    path = changed_alignment(dataset, tmp_path, {"s3": {column: derived}})
    findings = audit(written[0], dataset, alignment_path=path)

    homoplasy = of(findings, "homoplasy", PROBLEM)
    assert [finding.node for finding in homoplasy] == ["A.1"]
    assert "s3" in homoplasy[0].detail and "B.1 (1)" in homoplasy[0].detail


def test_swapped_alleles_are_reported_as_polarity(written, dataset):
    outdir, _scores = written
    column, ancestral, derived = PLANTED["A.1"]
    sequence, position = dataset.contig_of(column)
    edit(outdir, SNPS, f"\t{sequence}\t{position}\t{derived}\t{ancestral}\n",
         f"\t{sequence}\t{position}\t{ancestral}\t{derived}\n")
    # The isolates below A.1 now carry what the scheme calls ancestral, and
    # every genome outside it what the scheme calls derived.
    findings = audit(outdir, dataset)
    polarity = of(findings, "polarity", PROBLEM)
    assert [finding.node for finding in polarity] == ["A.1"]


def test_a_failed_snp_is_offered_a_replacement(written, dataset, tmp_path):
    # Spoil A.1's SNP in one isolate, and plant a second canonical site for A.1
    # elsewhere, which the audit should offer instead.
    column, ancestral, _derived = PLANTED["A.1"]
    records = read_alignment(dataset.alignment)
    spare = 650
    base = records["ref"][spare]
    other = next(letter for letter in "ACGT" if letter != base)
    changes = {"s1b": {column: ancestral}, "s1": {spare: other}}
    changes["s1b"][spare] = other
    path = changed_alignment(dataset, tmp_path, changes)

    replacement = of(audit(written[0], dataset, alignment_path=path), "replacement")
    assert [finding.node for finding in replacement] == ["A.1"]
    assert f"{spare + 1} (c1:{spare + 1}, {base}>{other})" in replacement[0].detail


def test_a_site_no_genome_has_a_base_at_is_reported(written, dataset, tmp_path):
    column, _ancestral, _derived = PLANTED["B.1"]
    path = changed_alignment(dataset, tmp_path, {"s4": {column: "N"}})
    missing = of(audit(written[0], dataset, alignment_path=path), "missing-data")
    assert [finding.node for finding in missing] == ["B.1"]
    assert "s4" in missing[0].detail


# -- more than one canSNP per node -------------------------------------------


def write_snp_set(outdir, extra_rows):
    """A canSNP set whose rank 1 is the SNP file, with these further rows."""
    with open(os.path.join(outdir, SNPS)) as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle][1:]
    lines = ["node\trank\tgenome\tsequence\tsequence_position\tderived_base\tancestral_base\t"
             "score"]
    lines += [f"{row[0]}\t1\t{row[3]}\t{row[4]}\t{row[5]}\t{row[6]}\t{row[7]}\t500"
              for row in rows]
    lines += extra_rows
    path = os.path.join(outdir, SNPS.replace("_snp.txt", "_canSNPs.txt"))
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    return rows[0][3]


def test_every_ranked_canSNP_is_checked_in_the_genomes(written, dataset):
    outdir = written[0]
    column, ancestral, derived, _isolates = INCONGRUENT
    genome = write_snp_set(outdir, [])
    sequence, position = dataset.contig_of(column)
    write_snp_set(outdir, [f"A.1\t2\t{genome}\t{sequence}\t{position}\t{derived}\t{ancestral}"
                           f"\t400"])
    scheme = read_scheme(outdir)
    assert scheme.problems == []
    assert [snp.rank for snp in scheme.checked_snps if snp.node == "A.1"] == [1, 2]

    findings = audit(outdir, dataset)
    # s1b, below A.1, lacks the derived base, and s3, outside it, carries it.
    for category in ("not-canonical", "homoplasy"):
        found = of(findings, category, PROBLEM)
        assert [finding.node for finding in found] == ["A.1"]
        assert found[0].detail.startswith("the canSNP of rank 2: ")
    # The SNP CanSNPer2 calls is sound, so no replacement is offered for it.
    assert of(findings, "replacement") == []


def test_a_rank_one_that_is_not_the_snp_file_is_a_problem(written, dataset):
    outdir = written[0]
    genome = write_snp_set(outdir, [])
    path = os.path.join(outdir, SNPS.replace("_snp.txt", "_canSNPs.txt"))
    with open(path) as handle:
        text = handle.read()
    column, _ancestral, _derived = PLANTED["B.1"]
    sequence, position = dataset.contig_of(column)
    text = text.replace(f"B.1\t1\t{genome}\t{sequence}\t{position}\t",
                        f"B.1\t1\t{genome}\t{sequence}\t{position + 1}\t")
    with open(path, "w") as handle:
        handle.write(text)
    problems = read_scheme(outdir).problems
    assert any("node 'B.1' a rank 1 that is not the SNP" in problem for problem in problems)


# -- extending the scheme ----------------------------------------------------


def test_genomes_the_scheme_does_not_resolve_are_a_new_clade(written, dataset, tmp_path):
    # Place s1 and s1b at A rather than at A.1. The site that marks them is
    # then a clade the scheme has no node for, below A.
    placements = {node: list(genomes) for node, genomes in CLUSTERS.items()}
    placements["A"] = placements.pop("A.1")
    groups = write_groups(tmp_path, placements)

    findings = audit(written[0], dataset, group_path=groups)
    clades = of(findings, "new-clade")
    assert [finding.node for finding in clades] == ["A"]
    assert clades[0].action == "extend"
    assert "s1, s1b" in clades[0].detail
    assert "A.1" in of(findings, "untested")[0].node


def test_a_clade_needs_enough_genomes_and_snps_to_be_reported(written, dataset, tmp_path):
    placements = {node: list(genomes) for node, genomes in CLUSTERS.items()}
    placements["A"] = placements.pop("A.1")
    groups = write_groups(tmp_path, placements)
    assert of(audit(written[0], dataset, group_path=groups, min_clade_size=3), "new-clade") == []
    assert of(audit(written[0], dataset, group_path=groups, min_clade_snps=2), "new-clade") == []


# -- the inputs RepGenR2 writes ----------------------------------------------


def test_genomes_named_by_file_name_are_matched_to_records(written, dataset, tmp_path):
    genome_map = tmp_path / "genomes_map.tsv"
    genome_map.write_text("".join(
        f"{genome}.fna.gz\t{node}\n" for node, genomes in CLUSTERS.items()
        for genome in genomes))
    findings = audit(written[0], dataset, group_path=str(genome_map))
    assert [finding for finding in findings if finding.severity == PROBLEM] == []


def test_a_projected_alignment_is_read_as_reference_coordinates(written, dataset, tmp_path):
    # x2fa, snippy and RepGenR2 project every row onto the reference, writing a
    # gap in the reference row where the reference was not aligned and
    # stopping where the last aligned block ends.
    records = read_alignment(dataset.alignment)
    projected = {}
    for name, sequence in records.items():
        row = sequence[:980]
        projected[name] = row[:20] + "-" * 10 + row[30:]
    path = str(tmp_path / "projected.fasta")
    write_fasta(path, list(projected.items()))

    findings = audit(written[0], dataset, alignment_path=path)
    assert [finding for finding in findings if finding.severity == PROBLEM] == []
    assert len(of(findings, "single-site")) == len(PLANTED) - 1

    snp_file = str(tmp_path / "snps.vcf")
    assert cli_run("extract", "--genome-map", dataset.genome_map, "--hierarchy", dataset.hierarchy,
                   "--alignment", path, "--reference", dataset.reference,
                   "-o", snp_file) == 0
    found = {record.node: (record.chrom, record.pos) for record in vcf_records(snp_file)
             if record.node}
    assert found["B.2"] == ("c2", 1)
    assert found["A.1"] == ("c1", PLANTED["A.1"][0] + 1)


# -- the command line --------------------------------------------------------


def test_the_genomes_are_checked_from_the_command_line(written, dataset, tmp_path):
    outdir, _scores = written
    report = tmp_path / "findings.tsv"
    assert cli_run("audit", "--scheme-dir", outdir, "--reference", dataset.reference,
                   "--alignment", dataset.alignment, "--genome-map", dataset.genome_map,
                   "-o", str(report)) == 0
    assert "single-site" in report.read_text()


def test_an_alignment_needs_a_group_file_and_a_reference(written, dataset):
    outdir, _scores = written
    assert audit_main(outdir, alignment=dataset.alignment) == 1
    assert audit_main(outdir, alignment=dataset.alignment,
                      reference=dataset.reference) == 1
    assert audit_main(outdir, genome_map=dataset.genome_map,
                      reference=dataset.reference) == 1


def test_a_recommendation_alone_exits_zero(written, dataset):
    outdir, _scores = written
    # The planted data leaves only recommendations: single sites, the
    # incongruent site, the ungrouped reference and the shallow nodes.
    assert audit_main(outdir, reference=dataset.reference, alignment=dataset.alignment,
                      genome_map=dataset.genome_map) == 0
    assert audit_main(outdir, reference=dataset.reference, alignment=dataset.alignment,
                      genome_map=dataset.genome_map, strict=True) == 2


# -- without an outgroup -----------------------------------------------------


def small_scheme(snps):
    """A root with three children, A holding three genomes and X and Y one each."""
    from cansnpscorer.scheme import Scheme, SchemeSnp, hierarchy_from_paths

    paths = [["root"], ["root", "A"], ["root", "X"], ["root", "Y"]]
    hierarchy, _conflicts, _problems = hierarchy_from_paths(paths)
    rows = [SchemeSnp(node=node, strain="t", reference="t", genome="REF", sequence="c1", position=position,
                      derived=derived, ancestral=ancestral)
            for node, position, derived, ancestral in snps]
    return Scheme(organism="small", paths=paths, snps=rows, hierarchy=hierarchy)


SMALL_CLUSTERS = {"A": ["a1", "a2", "a3"], "X": ["x1"], "Y": ["y1"]}


def test_a_node_holding_more_genomes_than_its_sisters_keeps_its_site():
    # Every genome is placed, so nothing is outside the root. The site at
    # column 2 separates A's three genomes from x1 and y1: that is A's canSNP,
    # not a new clade of the two fewer genomes.
    from cansnpscorer.alignment_audit import check_against_alignment
    from cansnpscorer.model import Reference

    reference = Reference.from_records([("c1", "ACGTACGT")])
    # x1 is the reference genome, so it comes first.
    alignment = {
        "x1": "ACGTACGT", "y1": "ACGTACGT",
        "a1": "ACTTACGT", "a2": "ACTTACGT", "a3": "ACTTACGT",
    }
    scheme = small_scheme([("A", 3, "T", "G")])

    findings = check_against_alignment(scheme, alignment, reference, SMALL_CLUSTERS)
    assert [finding for finding in findings if finding.category == "alignment"] == []
    assert of(findings, "new-clade") == []
    assert of(findings, "not-canonical") == [] and of(findings, "homoplasy") == []
    assert [finding.node for finding in of(findings, "single-site")] == ["A"], (
        "the site is credited to A")


def test_a_snp_every_genome_carries_is_reported_as_untested():
    from cansnpscorer.alignment_audit import check_against_alignment
    from cansnpscorer.model import Reference

    reference = Reference.from_records([("c1", "ACGTACGT")])
    alignment = dict.fromkeys(("a1", "a2", "a3", "x1", "y1"), "ACGTACGT")
    scheme = small_scheme([("root", 1, "A", "C")])

    findings = check_against_alignment(scheme, alignment, reference, SMALL_CLUSTERS)
    assert [finding.node for finding in of(findings, "no-outgroup")] == ["root"]


def test_a_grouping_against_a_better_supported_one_is_not_a_new_clade():
    # Three sites group a1 with a2, one site groups a2 with a3. No tree holds
    # both, so only the first is a clade to add; the second is counted as a
    # conflicting site.
    from cansnpscorer.alignment_audit import check_against_alignment
    from cansnpscorer.model import Reference

    reference = Reference.from_records([("c1", "AAAAAAAA")])
    alignment = {
        "x1": "AAAAAAAA", "y1": "AAAAAAAA",
        "a1": "CCCAAAAA", "a2": "CCCCAAAA", "a3": "AAACAAAA",
    }
    scheme = small_scheme([])

    findings = check_against_alignment(scheme, alignment, reference, SMALL_CLUSTERS)
    clades = of(findings, "new-clade")
    assert len(clades) == 1
    assert "a1, a2" in clades[0].detail and "3 SNPs" in clades[0].detail
    conflicting = of(findings, "conflicting-sites")
    assert len(conflicting) == 1
    assert "1 of them because they group genomes against a better supported" in (
        conflicting[0].detail)


# -- each genome's own path through the scheme -------------------------------


def test_a_sound_scheme_has_no_misplaced_or_discordant_genome(written, dataset):
    findings = audit(written[0], dataset)
    assert of(findings, "misplaced") == []
    assert of(findings, "discordant-genome") == []


def test_a_genome_placed_on_the_wrong_branch_is_reported_where_it_would_be_typed(
        written, dataset, tmp_path):
    # s2 carries the canSNPs of A.2 but the genome map puts it at A.1, as a
    # cluster spanning the two would.
    groups = write_groups(tmp_path, {"A.1": ["s1", "s1b", "s2"], "B.1": ["s3"], "B.2": ["s4"]})
    misplaced = of(audit(written[0], dataset, group_path=groups), "misplaced")
    assert [finding.node for finding in misplaced] == ["A.1"]
    assert "s2" in misplaced[0].detail and "typed at A.2" in misplaced[0].detail
    assert "at least 1 canSNP of A.2" in misplaced[0].detail


def test_a_genome_that_contradicts_the_scheme_at_many_nodes_is_reported(
        written, dataset, tmp_path):
    changes = {"s4": {PLANTED[node][0]: PLANTED[node][2] for node in ("A", "A.1", "A.2")}}
    path = changed_alignment(dataset, tmp_path, changes)
    findings = audit(written[0], dataset, alignment_path=path)
    discordant = of(findings, "discordant-genome")
    assert len(discordant) == 1 and "s4 (3)" in discordant[0].detail
    # It is not also reported as misplaced.
    assert of(findings, "misplaced") == []
    assert of(audit(written[0], dataset, alignment_path=path, genome_conflicts=4),
              "discordant-genome") == []


def test_one_homoplastic_site_in_one_genome_is_left_to_the_snp_checks(
        written, dataset, tmp_path):
    column, _ancestral, derived = PLANTED["A.1"]
    path = changed_alignment(dataset, tmp_path, {"s4": {column: derived}})
    findings = audit(written[0], dataset, alignment_path=path)
    assert of(findings, "discordant-genome") == [] and of(findings, "misplaced") == []
    assert [finding.node for finding in of(findings, "homoplasy")] == ["A.1"]


def test_the_cli_takes_the_genome_conflicts_threshold(written, dataset, tmp_path):
    changes = {"s4": {PLANTED[node][0]: PLANTED[node][2] for node in ("A", "A.1")}}
    path = changed_alignment(dataset, tmp_path, changes)
    output = tmp_path / "findings.tsv"
    cli_run("audit", "--scheme-dir", written[0], "--reference", dataset.reference,
            "--alignment", path, "--genome-map", dataset.genome_map,
            "--genome-conflicts", "2", "-o", str(output))
    assert "discordant-genome" in output.read_text()


def test_a_genome_at_a_node_the_scheme_lacks_tests_no_snp(written, dataset, tmp_path):
    # As after export --prune: the genome map still names a node the scheme
    # left out. Its genome carries the derived base of the nodes above it, and
    # is not evidence that their SNPs are homoplastic.
    outdir, _scores = written
    placements = {node: list(genomes) for node, genomes in CLUSTERS.items() if node != "B.2"}
    placements["B.2.1"] = ["s4"]
    findings = audit(outdir, dataset, group_path=write_groups(tmp_path, placements))
    assert not of(findings, "homoplasy") and not of(findings, "not-canonical")
    assert not of(findings, "missing-genome")
    unknown = of(findings, "unknown-cluster")
    assert len(unknown) == 1 and "1 genomes at 1 nodes" in unknown[0].detail
    assert all("s4" not in finding.detail for finding in of(findings, "ungrouped"))


# -- reading only the scheme's columns ---------------------------------------


def small_alignment(tmp_path, records=None):
    path = str(tmp_path / "small.fasta")
    write_fasta(path, records or [("ref", "ACGTACGT"), ("g1", "ACGAACGT"), ("g2", "TCGTACGA")])
    return path


def test_column_rows_match_full_read(tmp_path):
    from cansnpscorer.helpers.file_utils import read_alignment_columns

    path = small_alignment(tmp_path)
    full = read_alignment(path)
    kept = read_alignment_columns(path, lambda reference_row: {0, 4})
    assert list(kept) == list(full)
    assert kept["ref"] == full["ref"]
    assert all(kept[name][4] == full[name][4] and kept[name][0] == full[name][0]
               for name in full)
    assert all(len(kept[name]) == len(full[name]) for name in full)


def test_columns_are_asked_for_once_from_the_reference_row(tmp_path):
    from cansnpscorer.helpers.file_utils import read_alignment_columns

    path = small_alignment(tmp_path)
    seen = []
    read_alignment_columns(path, lambda row: seen.append(row) or {1})
    assert seen == ["ACGTACGT"]


def test_unread_column_raises(tmp_path):
    import pytest

    from cansnpscorer.helpers.file_utils import read_alignment_columns

    path = small_alignment(tmp_path)
    with pytest.raises(KeyError, match="column 2 was not read"):
        read_alignment_columns(path, lambda row: {0})["g1"][2]


def test_ragged_alignment_message_unchanged(tmp_path):
    import pytest

    from cansnpscorer.helpers.file_utils import read_alignment_columns

    path = small_alignment(tmp_path, [("ref", "ACGT"), ("g1", "ACG"), ("g2", "ACGT")])
    with pytest.raises(ValueError, match="the records are not aligned"):
        read_alignment_columns(path, lambda row: {0})
    with pytest.raises(ValueError, match="the records are not aligned"):
        read_alignment(path)


def test_a_single_record_and_a_repeated_name_are_refused(tmp_path):
    import pytest

    from cansnpscorer.helpers.file_utils import read_alignment_columns

    alone = small_alignment(tmp_path, [("ref", "ACGT")])
    with pytest.raises(ValueError, match="needs the reference and at least one other"):
        read_alignment_columns(alone, lambda row: {0})
    twice = small_alignment(tmp_path, [("ref", "ACGT"), ("g1", "ACGT"), ("g1", "ACGA")])
    with pytest.raises(ValueError, match="more than once"):
        read_alignment_columns(twice, lambda row: {0})
