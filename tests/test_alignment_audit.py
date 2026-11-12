"""The audit against the genomes themselves.

The scheme is the one the export subcommand writes from the planted dataset, so
every SNP in it is canonical in that dataset's alignment. Each test changes the
alignment, the group file or the scheme in one way and asserts what the audit
says about it.
"""

from conftest import CLUSTERS, PLANTED, SNPS, cli_run, edit, write_fasta

from cansnpscorer.cansnp_audit import PROBLEM, audit_scheme
from cansnpscorer.cansnp_audit import main as audit_main
from cansnpscorer.helpers.file_utils import (
    HCParser,
    read_alignment,
    recognize_group_file_format,
    reference_parser,
)
from cansnpscorer.scheme import read_scheme


def audit(outdir, dataset, alignment_path=None, group_path=None, **arguments):
    alignment = read_alignment(alignment_path or dataset.alignment)
    group_path = group_path or dataset.groups_csv
    clusters = HCParser(group_path, recognize_group_file_format(group_path)).isolates
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
    path = tmp_path / "groups.csv"
    path.write_text("Strain,Group\n" + "".join(
        f"{genome},{node}\n" for node, genomes in placements.items() for genome in genomes))
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
    edit(outdir, SNPS, f"\t{column + 1}\t{derived}\t{ancestral}\n",
         f"\t{column + 1}\t{ancestral}\t{derived}\n")
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


def test_a_repgenr2_genome_map_places_the_genomes(written, dataset, tmp_path):
    genome_map = tmp_path / "genomes_map.tsv"
    genome_map.write_text("".join(
        f"{genome}.fna.gz\t{node}\n" for node, genomes in CLUSTERS.items()
        for genome in genomes))
    assert recognize_group_file_format(str(genome_map)) == "genome_map"
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

    snp_file = str(tmp_path / "snps.tsv")
    assert cli_run("extract", "--groupfile", dataset.groups_csv, "--hcfile", dataset.hc_yaml,
                   "--fastafile", path, "--reference-file", dataset.reference,
                   "-o", snp_file) == 0
    with open(snp_file) as handle:
        rows = [line.split("\t") for line in handle.read().splitlines()[1:]]
    found = {row[2]: (row[0], int(row[1])) for row in rows if row[2] != "None"}
    assert found["B.2"] == ("c2", 1)
    assert found["A.1"] == ("c1", PLANTED["A.1"][0] + 1)


# -- the command line --------------------------------------------------------


def test_the_genomes_are_checked_from_the_command_line(written, dataset, tmp_path):
    outdir, _scores = written
    report = tmp_path / "findings.tsv"
    assert cli_run("audit", "--scheme-dir", outdir, "--reference-file", dataset.reference,
                   "--fastafile", dataset.alignment, "--groupfile", dataset.groups_csv,
                   "-o", str(report)) == 0
    assert "single-site" in report.read_text()


def test_an_alignment_needs_a_group_file_and_a_reference(written, dataset):
    outdir, _scores = written
    assert audit_main(outdir, fastafile=dataset.alignment) == 1
    assert audit_main(outdir, fastafile=dataset.alignment,
                      reference_file=dataset.reference) == 1
    assert audit_main(outdir, groupfile=dataset.groups_csv,
                      reference_file=dataset.reference) == 1


def test_a_recommendation_alone_exits_zero(written, dataset):
    outdir, _scores = written
    # The planted data leaves only recommendations: single sites, the
    # incongruent site, the ungrouped reference and the shallow nodes.
    assert audit_main(outdir, reference_file=dataset.reference, fastafile=dataset.alignment,
                      groupfile=dataset.groups_csv) == 0
    assert audit_main(outdir, reference_file=dataset.reference, fastafile=dataset.alignment,
                      groupfile=dataset.groups_csv, strict=True) == 2


# -- without an outgroup -----------------------------------------------------


def small_scheme(snps):
    """A root with three children, A holding three genomes and X and Y one each."""
    from cansnpscorer.scheme import Scheme, SchemeSnp, hierarchy_from_paths

    paths = [["root"], ["root", "A"], ["root", "X"], ["root", "Y"]]
    hierarchy, _conflicts, _problems = hierarchy_from_paths(paths)
    rows = [SchemeSnp(node=node, strain="t", reference="t", genome="REF", position=position,
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
    alignment = {name: "ACGTACGT" for name in ("a1", "a2", "a3", "x1", "y1")}
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
