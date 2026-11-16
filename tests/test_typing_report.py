"""Typing the genomes of an alignment with a scheme, and what is made of it:
the typing report of audit, the ranking by observed behaviour, and the choice
of a node's canSNPs so that they fail independently."""

import csv
import os

from conftest import CLUSTERS, PLANTED, SNPS, cli_run, vcf_records, write_fasta

from cansnpscorer.helpers.data_utils import ranking, select_per_node
from cansnpscorer.helpers.file_utils import read_alignment
from cansnpscorer.helpers.scoring_methods import ScoredSnp
from cansnpscorer.typing_report import (
    ANCESTRAL,
    DERIVED,
    MISSING,
    node_state,
    read_observed,
    set_call,
    single_call,
)

LINEAGE = {"root": ["root"], "A": ["A", "root"], "A.1": ["A.1", "A", "root"],
           "A.2": ["A.2", "A", "root"], "B": ["B", "root"]}


def table(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_report(written, dataset, tmp_path, alignment=None, *extra):
    outdir, _scores = written
    prefix = str(tmp_path / "typing")
    status = cli_run("audit", "--scheme-dir", outdir, "--reference", dataset.reference,
                     "--alignment", alignment or dataset.alignment,
                     "--typing-report", prefix, *extra)
    return status, prefix


def changed(dataset, tmp_path, changes):
    records = read_alignment(dataset.alignment)
    for name, columns in changes.items():
        sequence = list(records[name])
        for column, base in columns.items():
            sequence[column] = base
        records[name] = "".join(sequence)
    path = str(tmp_path / "changed.fasta")
    write_fasta(path, list(records.items()))
    return path


def add_set(outdir, extra):
    """Write the scheme's canSNP set: its SNP file as rank 1, and further
    ranks given as {node: [(column, derived, ancestral), ...]}."""
    with open(os.path.join(outdir, SNPS)) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    lines = ["node\trank\tgenome\tposition\tderived_base\tancestral_base\tscore"]
    for row in rows:
        lines.append("\t".join([row["snp_id"], "1", row["genome"], row["position"],
                                row["derived_base"], row["ancestral_base"], "600"]))
        for rank, (column, derived, ancestral) in enumerate(extra.get(row["snp_id"], []), 2):
            lines.append("\t".join([row["snp_id"], str(rank), row["genome"], str(column + 1),
                                    derived, ancestral, "600"]))
    with open(os.path.join(outdir, SNPS.replace("_snp.txt", "_canSNPs.txt")), "w") as handle:
        handle.write("\n".join(lines) + "\n")


# -- the two rules -----------------------------------------------------------


def test_a_node_takes_the_state_of_most_of_its_cansnps():
    assert node_state([DERIVED, DERIVED, ANCESTRAL]) == DERIVED
    assert node_state([DERIVED, ANCESTRAL, ANCESTRAL]) == ANCESTRAL
    assert node_state([DERIVED, ANCESTRAL]) == MISSING
    assert node_state([MISSING]) == MISSING


def test_the_single_rule_needs_a_complete_chain_of_enough_nodes():
    states = {"root": DERIVED, "A": DERIVED, "A.1": DERIVED}
    assert single_call(states, LINEAGE, 3) == "A.1"
    assert single_call({**states, "A": MISSING}, LINEAGE, 3) is None
    assert single_call({"root": DERIVED, "A": DERIVED}, LINEAGE, 3) is None
    assert single_call({"root": DERIVED, "A": DERIVED}, LINEAGE, 2) == "A"


def test_the_set_rule_passes_a_node_of_unknown_state_but_not_an_ancestral_one():
    assert set_call({"root": DERIVED, "A": MISSING, "A.1": DERIVED}, LINEAGE) == "A.1"
    assert set_call({"root": DERIVED, "A": ANCESTRAL, "A.1": DERIVED}, LINEAGE) == "root"


def test_two_deepest_nodes_on_different_branches_give_no_call():
    states = {"root": DERIVED, "A": DERIVED, "A.1": DERIVED, "A.2": DERIVED}
    assert set_call(states, LINEAGE) is None
    assert single_call(states, LINEAGE, 3) is None


# -- the report --------------------------------------------------------------


def test_every_isolate_is_typed_at_its_node_under_both_rules(written, dataset, tmp_path):
    status, prefix = run_report(written, dataset, tmp_path, None,
                                "--genome-map", dataset.genome_map)
    assert status in (0, 2)  # the audit's own findings do not matter here
    genomes = {row["genome"]: row for row in table(prefix + "_genomes.tsv")}
    for node, isolates in CLUSTERS.items():
        for isolate in isolates:
            row = genomes[isolate]
            assert (row["single_call"], row["set_call"]) == (node, node)
            assert (row["single_outcome"], row["set_outcome"]) == ("same", "same")
    assert genomes["ref"]["set_call"] == "NA"

    nodes = {row["node"]: row for row in table(prefix + "_nodes.tsv")}
    assert nodes["A"]["through"] == "3" and nodes["A.1"]["set_calls"] == "2"
    summary = {row["metric"]: row for row in table(prefix + "_summary.tsv")}
    assert summary["typed"]["single"] == summary["typed"]["set"] == "5"
    assert summary["call is the expected node"]["set"] == "5"


def test_the_report_needs_no_placement_of_the_genomes(written, dataset, tmp_path):
    status, prefix = run_report(written, dataset, tmp_path)
    assert status in (0, 2)
    rows = table(prefix + "_genomes.tsv")
    assert {row["expected"] for row in rows} == {"NA"}
    assert sum(row["set_call"] != "NA" for row in rows) == 5


def test_without_the_report_an_alignment_still_needs_a_placement(written, dataset):
    outdir, _scores = written
    assert cli_run("audit", "--scheme-dir", outdir, "--reference", dataset.reference,
                   "--alignment", dataset.alignment) == 1


def test_a_missing_base_stops_the_single_rule_and_not_the_set_rule(written, dataset,
                                                                   tmp_path):
    alignment = changed(dataset, tmp_path, {"s1": {PLANTED["A"][0]: "N"}})
    _status, prefix = run_report(written, dataset, tmp_path, alignment)
    row = {r["genome"]: r for r in table(prefix + "_genomes.tsv")}["s1"]
    assert row["single_call"] == "NA"
    assert (row["single_stop_node"], row["single_stop_reason"]) == ("A", "missing")
    assert row["set_call"] == "A.1" and row["missing_on_path"] == "1"
    nodes = {r["node"]: r for r in table(prefix + "_nodes.tsv")}
    assert nodes["A"]["single_stops"] == "1" and nodes["A"]["unknown_on_path"] == "1"
    snps = {r["node"]: r for r in table(prefix + "_snps.tsv")}
    assert snps["A"]["missing"] == "1" and snps["A"]["called"] == "5"


def test_a_set_outvotes_one_discordant_base_and_the_snp_is_counted(written, dataset,
                                                                    tmp_path):
    outdir, _scores = written
    column, ancestral, _derived = PLANTED["A"]
    # Two more canSNPs for A, whose ancestral base is the reference's own.
    reference = read_alignment(dataset.alignment)["ref"]
    added = [(at, "C" if reference[at] != "C" else "G", reference[at]) for at in (210, 220)]
    extra = {name: {at: new if name in ("s1", "s1b", "s2") else old
                    for at, new, old in added}
             for name in read_alignment(dataset.alignment)}
    extra["s2"][column] = ancestral  # the first canSNP of A fails in s2
    alignment = changed(dataset, tmp_path, extra)
    add_set(outdir, {"A": added})

    _status, prefix = run_report(written, dataset, tmp_path, alignment)
    row = {r["genome"]: r for r in table(prefix + "_genomes.tsv")}["s2"]
    assert row["single_call"] == "NA"
    assert (row["single_stop_node"], row["single_stop_reason"]) == ("A", "ancestral")
    assert row["set_call"] == "A.2"

    snps = {(r["node"], r["rank"]): r for r in table(prefix + "_snps.tsv")}
    assert snps[("A", "1")]["discordant_in_set"] == "1"
    assert snps[("A", "1")]["ancestral_on_path"] == "1"
    assert snps[("A", "2")]["discordant_in_set"] == "0"
    nodes = {r["node"]: r for r in table(prefix + "_nodes.tsv")}
    assert nodes["A"]["snps"] == "3" and nodes["A"]["min_margin"] == "1"
    assert nodes["A"]["split"] == "1" and nodes["A"]["all_derived"] == "2"

    observed = read_observed(prefix + "_snps.tsv")
    assert observed[("A", column + 1)].discordant == 1
    assert observed[("A", 211)].discordant == 0


# -- ranking by what was observed -------------------------------------------


def snp(position, total=600, contig="c1", locus="", **fields):
    return ScoredSnp(contig=contig, position=position, node="A", total=total, terms={},
                     derived="C", ancestral="G", locus=locus, **fields)


def test_a_snp_that_disagreed_with_its_node_is_ranked_after_a_lower_score():
    failed, sound = snp(10, total=600, discordant=1), snp(20, total=300)
    assert ranking(sound) < ranking(failed)
    assert select_per_node([failed, sound], top=1)["A"] == [sound]


def test_observed_missing_counts_with_the_missing_of_the_alignment():
    assert ranking(snp(10, total=300)) < ranking(snp(20, observed_missing=1))
    assert ranking(snp(10, missing=1)) < ranking(snp(20, missing=1, observed_missing=1))


def test_evaluate_reads_the_observed_behaviour_into_the_scores(written, dataset, tmp_path):
    _status, prefix = run_report(written, dataset, tmp_path,
                                 changed(dataset, tmp_path, {"s1": {PLANTED["A"][0]: "N"}}))
    snp_file, scores = str(tmp_path / "snps.vcf"), str(tmp_path / "observed.vcf")
    assert cli_run("evaluate", "--snps", snp_file, "--reference", dataset.reference,
                   "--genes", dataset.gff, "--repeats", dataset.repeats, "-o", scores,
                   "--profile", "typing", "--observed", prefix + "_snps.tsv") == 0
    by_node = {record.node: record for record in vcf_records(scores)}
    assert by_node["A"].info["OBS_MISSING"] == "1"
    assert by_node["B"].info["OBS_MISSING"] == "0"
    assert by_node["B"].info["OBS_DISCORDANT"] == "0"


def test_the_typing_profile_is_the_default_and_drops_the_gc_term_and_the_cutoff(
        written, dataset, tmp_path):
    snp_file, scores = str(tmp_path / "snps.vcf"), str(tmp_path / "typing.vcf")
    assert cli_run("evaluate", "--snps", snp_file, "--reference", dataset.reference,
                   "--genes", dataset.gff, "--repeats", dataset.repeats, "-o", scores) == 0
    records = vcf_records(scores)
    assert records and all(record.info["GC"] == "0" for record in records)
    assert all("LowScore" not in record.filters for record in records)
    assert all(int(record.info["GENE"]) <= 100 for record in records)


def test_the_assay_profile_keeps_its_cutoff_of_500(written, dataset, tmp_path):
    snp_file, scores = str(tmp_path / "snps.vcf"), str(tmp_path / "assay.vcf")
    assert cli_run("evaluate", "--snps", snp_file, "--reference", dataset.reference,
                   "--genes", dataset.gff, "--repeats", dataset.repeats, "-o", scores,
                   "--profile", "assay") == 0
    records = vcf_records(scores)
    assert any(int(record.info["GC"]) > 0 for record in records)
    assert all(("LowScore" in record.filters) == (int(record.info["SCORE"]) <= 500)
               for record in records)


# -- sets that fail independently -------------------------------------------


def test_the_first_of_a_set_is_the_best_ranked_whatever_its_locus():
    best, other = snp(10, total=600, locus="1-900"), snp(5000, total=500, locus="4000-6000")
    assert select_per_node([other, best], top=2)["A"] == [best, other]


def test_another_gene_is_taken_before_a_better_snp_of_the_same_gene():
    first = snp(100, total=600, locus="1-9000")
    same_gene = snp(5000, total=590, locus="1-9000")
    other_gene = snp(20000, total=400, locus="19000-21000")
    assert select_per_node([first, same_gene, other_gene], top=2)["A"] == [first, other_gene]
    assert select_per_node([first, same_gene, other_gene], top=3)["A"] == [
        first, other_gene, same_gene]


def test_of_equal_scores_an_unused_contig_is_taken_first():
    first = snp(100, locus="1-900")
    same_contig = snp(5000, locus="4000-6000")
    other_contig = snp(7000, contig="c2", locus="6500-7500")
    assert select_per_node([first, same_contig, other_contig], top=2)["A"] == [
        first, other_contig]
    lower = snp(7000, total=599, contig="c2", locus="6500-7500")  # score comes first
    assert select_per_node([first, same_contig, lower], top=2)["A"] == [first, same_contig]


def test_a_less_reliable_snp_is_not_taken_for_its_locus():
    first = snp(100, locus="1-9000")
    same_gene = snp(5000, locus="1-9000")
    missing = snp(20000, locus="19000-21000", missing=1)
    assert select_per_node([first, same_gene, missing], top=2)["A"] == [first, same_gene]


def test_spacing_still_holds_between_genes_and_outside_them():
    first, near = snp(1000, locus="1-1200"), snp(1500, locus="1300-2000")
    far = snp(9000, total=100)
    assert select_per_node([first, near, far], top=3, spacing=1000)["A"] == [first, far]
    assert select_per_node([first, near, far], spacing=0)["A"] == [first, near, far]


def test_the_locus_of_a_snp_survives_the_scores_file(written, dataset):
    from cansnpscorer import vcf
    from cansnpscorer.helpers.file_utils import reference_parser

    _outdir, scores = written
    read = vcf.read_scores(scores, reference_parser(dataset.reference))
    in_gene = [score for score in read if score.terms["gene_score"] > 0]
    assert in_gene and all(score.locus for score in in_gene)
    assert all(not score.locus for score in read if score.terms["gene_score"] == 0)


def test_a_snp_that_was_typed_with_and_held_is_ranked_before_one_never_observed():
    held, untried = snp(10, total=300, observed_genomes=50), snp(20, total=600)
    failed = snp(30, total=600, observed_genomes=50, discordant=2)
    assert sorted([failed, untried, held], key=ranking) == [held, untried, failed]


def test_what_was_observed_survives_the_scores_file(written, dataset, tmp_path):
    from cansnpscorer import vcf
    from cansnpscorer.helpers.file_utils import reference_parser

    _status, prefix = run_report(written, dataset, tmp_path,
                                 changed(dataset, tmp_path, {"s1": {PLANTED["A"][0]: "N"}}))
    scores = str(tmp_path / "observed.vcf")
    assert cli_run("evaluate", "--snps", str(tmp_path / "snps.vcf"), "--reference",
                   dataset.reference, "--genes", dataset.gff, "--repeats", dataset.repeats,
                   "-o", scores, "--profile", "typing",
                   "--observed", prefix + "_snps.tsv") == 0
    read = {score.node: score for score in
            vcf.read_scores(scores, reference_parser(dataset.reference))}
    assert (read["A"].observed_genomes, read["A"].observed_missing) == (5, 1)
    assert (read["B"].observed_genomes, read["B"].observed_missing) == (6, 0)


def test_audit_findings_and_typing_unchanged_when_only_scheme_columns_are_read(
        written, dataset, tmp_path):
    outdir, _scores = written
    outputs = {}
    for label, extra in (("full", []), ("columns", ["--no-new-clades"])):
        prefix = str(tmp_path / label)
        findings = str(tmp_path / f"{label}_findings.tsv")
        status = cli_run("audit", "--scheme-dir", outdir, "--reference", dataset.reference,
                         "--alignment", dataset.alignment, "--genome-map", dataset.genome_map,
                         "--typing-report", prefix, "-o", findings, *extra)
        outputs[label] = (status, open(findings).read(),
                          open(prefix + "_genomes.tsv").read(),
                          open(prefix + "_snps.tsv").read())
    assert outputs["full"][0] == outputs["columns"][0]
    # Without the search for new clades, the findings are the others of the full run.
    assert set(outputs["columns"][1].splitlines()) <= set(outputs["full"][1].splitlines())
    assert outputs["full"][2:] == outputs["columns"][2:]


def test_typing_only_audit_reads_the_scheme_columns(written, dataset, tmp_path, caplog):
    import logging

    caplog.set_level(logging.INFO)
    status, prefix = run_report(written, dataset, tmp_path)
    assert status in (0, 2)
    assert "scheme's columns" in caplog.text


# -- groups of a node's canSNPs by the genomes they are derived in -----------


def derived_set_rows(states_by_genome, discordant=()):
    """The _snps.tsv rows of one node 'A' whose canSNPs have these states,
    given per genome as one character per canSNP."""
    from cansnpscorer.model import Reference
    from cansnpscorer.scheme import SchemeSnp
    from cansnpscorer.typing_report import GenomeTyping, _snp_rows

    width = len(next(iter(states_by_genome.values())))
    snps = [(SchemeSnp("A", "", "", "G", 10 + index, "T", "C", rank=index + 1), index)
            for index in range(width)]
    typed = [GenomeTyping(genome=name, states=states, discordant=name in discordant)
             for name, states in states_by_genome.items()]
    reference = Reference.from_records([("c1", "A" * 400)])
    return _snp_rows({"A": snps}, {"A": 0}, typed, reference)


def labels_of(rows):
    from cansnpscorer.typing_report import SNP_COLUMNS

    at = SNP_COLUMNS.index("derived_set")
    return [int(row[at]) if row[at] else None for row in rows]


def test_derived_set_is_the_last_column_of_the_snp_table():
    from cansnpscorer.typing_report import SNP_COLUMNS

    assert SNP_COLUMNS[-1] == "derived_set"
    assert SNP_COLUMNS[:-1][-3:] == ("discordant_in_set", "ancestral_on_path",
                                     "derived_off_path")


def test_equal_counts_different_genomes_are_two_groups():
    # Two canSNPs derived in g1-g3, two in g4-g6: the counts are equal.
    rows = derived_set_rows({"g1": "DDAA", "g2": "DDAA", "g3": "DDAA",
                             "g4": "AADD", "g5": "AADD", "g6": "AADD", "g7": "AAAA"})
    assert labels_of(rows) == [1, 1, 2, 2]


def test_missing_base_does_not_split_a_group():
    # The first is derived in g1 and g2 and has no base in g3.
    rows = derived_set_rows({"g1": "DDD", "g2": "DDD", "g3": "NDD", "g4": "AAA"})
    assert labels_of(rows) == [1, 1, 1]


def test_the_set_of_most_genomes_is_labelled_first_and_none_without_derived_genomes():
    rows = derived_set_rows({"g1": "DAA", "g2": "ADA", "g3": "ADA", "g4": "ADA",
                             "g5": "ADA", "g6": "AAA"})
    assert labels_of(rows) == [2, 1, None]


def test_a_discordant_genome_does_not_split_a_group():
    rows = derived_set_rows({"g1": "DD", "g2": "DD", "bad": "DA", "g4": "AA"},
                            discordant={"bad"})
    assert labels_of(rows) == [1, 1]


def test_a_node_of_too_many_cansnps_is_labelled_for_the_first_ones_only(caplog):
    import logging

    from cansnpscorer.typing_report import MAX_SET_LABELLED

    width = MAX_SET_LABELLED + 3
    with caplog.at_level(logging.INFO):
        rows = derived_set_rows({"g1": "D" * width, "g2": "D" * width, "g3": "A" * width})
    assert labels_of(rows) == [1] * MAX_SET_LABELLED + [None] * 3
    assert "derived_set" in caplog.text and "A" in caplog.text


def group_table(tmp_path, rows, with_sets=True):
    """A _snps.tsv table of node 'A' from (derived, missing, derived_set) per canSNP,
    of 300 genomes."""
    from cansnpscorer.typing_report import SNP_COLUMNS

    columns = tuple(column for column in SNP_COLUMNS if with_sets or column != "derived_set")
    lines = ["\t".join(columns)]
    for index, (derived, missing, label) in enumerate(rows):
        fields = ["A", str(index + 1), "G", str(100 + index), "c1", str(100 + index), "T", "C",
                  "300", str(300 - missing), str(missing), str(derived),
                  str(300 - missing - derived), "0", "5", "0", "0"]
        if with_sets:
            fields.append(str(label) if label else "")
        lines.append("\t".join(fields))
    path = tmp_path / "group_snps.tsv"
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def test_derived_set_groups_canSNPs_of_equal_counts_apart(tmp_path):
    from cansnpscorer.typing_report import judge_by_narrowest_group, width_groups

    # Three sets of two canSNPs each, two of them of equal count: by counts
    # four of six canSNPs are one group, by sets none holds a majority.
    rows = [(150, 0, 1), (150, 0, 1), (60, 0, 2), (60, 0, 2), (60, 0, 3), (59, 0, 3)]
    observed = read_observed(group_table(tmp_path, rows))
    assert [row.derived_set for row in observed.values()] == [1, 1, 2, 2, 3, 3]
    assert [len(group) for group in width_groups(list(observed.values()))] == [2, 2, 2]
    judged, nodes = judge_by_narrowest_group(observed)
    assert nodes == ["A"]
    # The narrowest group is the last set; its most common count is 59 (ties
    # go to the smaller count).
    assert [judged[("A", 100 + index)].discordant for index in range(6)] == [91, 91, 1, 1, 1, 0]


def test_old_table_without_derived_set_reads_and_groups_by_count(tmp_path):
    from cansnpscorer.typing_report import judge_by_narrowest_group, width_groups

    rows = [(150, 0, 0), (150, 0, 0), (60, 0, 0), (60, 0, 0), (60, 0, 0), (59, 0, 0)]
    observed = read_observed(group_table(tmp_path, rows, with_sets=False))
    assert all(row.derived_set == 0 for row in observed.values())
    assert [len(group) for group in width_groups(list(observed.values()))] == [2, 4]
    judged, nodes = judge_by_narrowest_group(observed)
    assert nodes == [] and judged == observed


def test_a_node_with_an_unlabelled_derived_cansnp_groups_by_count(tmp_path):
    from cansnpscorer.typing_report import width_groups

    rows = [(150, 0, 1), (150, 0, 1), (60, 0, 2), (60, 0, 2), (60, 0, 0), (59, 0, 3)]
    observed = read_observed(group_table(tmp_path, rows))
    assert [len(group) for group in width_groups(list(observed.values()))] == [2, 4]


def test_derived_set_allows_conflicts_up_to_the_width_tolerance():
    # 100 genomes, tolerance 2. The second canSNP is ancestral in one genome
    # the first is derived in; the third in three.
    states = {f"g{number}": "DDD" if number < 40 else "AAA" for number in range(100)}
    states["g0"] = "DAD"
    for number in (1, 2):
        states[f"g{number}"] = "DDA"
    states["g40"] = "ADD"
    assert labels_of(derived_set_rows(states))[:2] == [1, 1]
    assert labels_of(derived_set_rows(states))[2] == 2


def test_derived_set_a_cansnp_without_enough_bases_founds_no_label():
    # 10 genomes, tolerance 2. The first canSNP is derived in g1-g6 and has
    # no base in g7-g10; it would found label 1 and keep the other two apart.
    states = {"g1": "DDD", "g2": "DDD", "g3": "DDD", "g4": "DAA", "g5": "DAA", "g6": "DAA",
              "g7": "NDD", "g8": "NDD", "g9": "NDD", "g10": "NAA"}
    assert labels_of(derived_set_rows(states)) == [2, 1, 1]
