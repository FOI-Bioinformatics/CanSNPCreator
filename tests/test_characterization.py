"""Behaviour that is expected to survive the hardening work.

These tests describe what the tools do today and are intended to keep doing:
the parsers, the five score terms and their scaling, and the output columns.
Tests for the defects being fixed live beside them, one file per area.
"""

import pytest
from conftest import CLUSTERS, GENOMES, INCONGRUENT, PLANTED, node_isolates

from cansnpscorer.helpers.data_utils import (
    canSNP_extractor,
    codon_comparer,
    codon_complement,
    filter_scores_by_cutoff,
    gc_content,
    neighbour_distance,
    search,
    snp_extractor,
    summarize_scores,
)
from cansnpscorer.helpers.file_utils import (
    HCParser,
    clean_name,
    fasta_parser,
    gff_parser,
    read_hc_file,
    recognize_group_file_format,
    recognize_hc_file_format,
    score_file_writer,
    snp_file_parser,
    trf_parser,
)
from cansnpscorer.helpers.genetic_code import codon_table
from cansnpscorer.helpers.scoring_methods import MAXIMUM_SCORE, TERMS, ScoredSnp

# -- parsers -----------------------------------------------------------------


def test_fasta_parser_output_formats(dataset):
    as_dict = fasta_parser(dataset.alignment, "dict")
    assert list(as_dict) == GENOMES
    assert as_dict["ref"] == dataset.sequences["ref"]

    as_named = fasta_parser(dataset.alignment, "named_list")
    assert [name for name, _ in as_named] == GENOMES

    as_list = fasta_parser(dataset.alignment, "list")
    assert as_list[0] == dataset.sequences["ref"]


def test_clean_name_keeps_only_alphanumeric():
    assert clean_name("GCA_000007125.1 ASM712v1") == "GCA0000071251ASM712v1"


def test_format_recognition(dataset):
    assert recognize_group_file_format(dataset.groups_csv) == "csv"
    assert recognize_group_file_format(dataset.groups_repgenr) == "repgenr_group"
    assert recognize_hc_file_format(dataset.hc_yaml) == "yaml"
    assert recognize_hc_file_format(dataset.hc_repgenr) == "repgenr"


def test_gff_parser_reports_genes_per_region_zero_based(dataset):
    genes = gff_parser(dataset.gff)
    assert sorted(genes) == [1, 2]
    assert genes[1] == [[150, 449, "+"], [550, 649, "-"]]
    assert genes[2] == [[10, 199, "+"]]


def test_trf_parser_reports_intervals_per_contig(dataset):
    repeats = trf_parser(dataset.trf)
    assert repeats == {1: [(60, 95)], 2: [(120, 160)]}


def test_group_files_describe_the_same_clusters(dataset):
    from_csv = HCParser(dataset.groups_csv, data_format="csv").isolates
    from_repgenr = HCParser(dataset.groups_repgenr, data_format="repgenr_group").isolates
    assert from_csv == CLUSTERS
    assert from_repgenr == CLUSTERS


# -- extraction --------------------------------------------------------------


def test_snp_extractor_finds_the_biallelic_columns(dataset):
    snps = snp_extractor(fasta_parser(dataset.alignment, "dict"))
    found = [snp["index"] for snp in snps]
    expected = sorted([column for column, _, _ in PLANTED.values()] + [INCONGRUENT[0]])
    assert found == expected, "the tri-allelic and gapped columns must be left out"


def test_cansnp_extractor_assigns_every_planted_node(dataset):
    alignment = fasta_parser(dataset.alignment, "dict")
    snps = snp_extractor(alignment)
    tree = read_hc_file(dataset.hc_yaml, "yaml")
    groups = HCParser(dataset.groups_csv, data_format="csv").reader(tree)

    cansnps = canSNP_extractor(alignment, groups, snps)
    assigned = {int(position): node for position, node, *_ in cansnps}

    for node, (column, _ancestral, _derived) in PLANTED.items():
        assert assigned[column] == node
    assert assigned[INCONGRUENT[0]] == "None", "a SNP crossing two branches is not canonical"


def test_hierarchy_from_yaml_gives_internal_nodes_their_isolates(dataset):
    tree = read_hc_file(dataset.hc_yaml, "yaml")
    groups = HCParser(dataset.groups_csv, data_format="csv").reader(tree)
    for node in ["root", "A", "B", "A.1", "A.2", "B.1", "B.2"]:
        assert sorted(groups[node]) == sorted(node_isolates(node)), node


# -- scoring -----------------------------------------------------------------


def test_neighbour_distance_normalises_to_the_closest_neighbour():
    assert neighbour_distance(100, 90, 150) == pytest.approx(10 / 25)
    assert neighbour_distance(100, 0, 1000) == pytest.approx(1.0)
    assert neighbour_distance(100, 100, 150) == pytest.approx(0.0)


def test_gc_content_scores_distance_from_the_goal():
    assert gc_content(100, "AGCT" * 100, 0.5) == pytest.approx(1.0), "at the goal"
    assert gc_content(100, "GC" * 200, 0.5) == pytest.approx(0.0), "all GC"
    assert gc_content(100, "AT" * 200, 0.5) == pytest.approx(0.0), "no GC"


def test_search_reports_inside_near_and_absent():
    # The interval comes back only for a position inside it, where its reading
    # frame applies; near and absent both report a distance alone.
    intervals = [(100, 200), (400, 500)]
    assert search(150, intervals, 50) == (0.0, (100, 200))
    assert search(80, intervals, 50) == (pytest.approx(0.4), None)
    assert search(300, intervals, 50) == (1.0, None)


def test_codon_helpers():
    assert len(codon_table) == 64
    assert codon_complement("ATG") == "TAC"
    assert codon_comparer(["TTA", "TTG"], "+") == 1.0, "synonymous"
    assert codon_comparer(["TTA", "TGA"], "+") == 0.0, "introduces a stop"
    assert codon_comparer(["TTA", "GTA"], "+") == 0.5, "changes the amino acid"


def test_score_terms_are_scaled_and_summed(scored):
    for entry in scored:
        assert isinstance(entry.position, int)
        assert isinstance(entry.node, str)
        assert set(entry.terms) == set(TERMS)
        assert entry.total == sum(entry.terms.values())
        assert all(0 <= term <= 200 for term in entry.terms.values())


def test_gene_term_is_weighted_twice_by_default(scored):
    gene_terms = {entry.terms["gene_score"] for entry in scored}
    assert gene_terms <= {0, 200}, "the gene term is doubled, so it is 0 or 200"
    assert MAXIMUM_SCORE == 600


# -- filtering and output ----------------------------------------------------


def make_score(position, node, total):
    return ScoredSnp(contig="c1", position=position, node=node, total=total,
                     terms=dict.fromkeys(TERMS, 0))


def test_filter_scores_by_cutoff():
    scores = [make_score(1, "A", 300), make_score(2, "A", 100), make_score(3, "B", 400)]
    assert [entry.position for entry in filter_scores_by_cutoff(scores, 150)] == [1, 3]
    top = filter_scores_by_cutoff(scores, 50, top=1)
    assert sorted(entry.position for entry in top) == [1, 3]


def test_summarize_scores_counts_snps_per_node():
    scores = [make_score(1, "A", 300), make_score(2, "A", 300), make_score(3, "B", 300)]
    assert summarize_scores(["A", "B", "C"], scores) == {"A": 2, "B": 1}


def test_score_file_writer_columns(tmp_path, scored):
    output = str(tmp_path / "scores.csv")
    score_file_writer(scored, output)
    lines = open(output).read().splitlines()
    assert lines[0] == (
        "Contig,Position,Group,Score,NeighborDistanceScore,GCContentScore,"
        "TandemRepeatProximityScore,GeneProximityScore,MutationTypeScore"
    )
    assert len(lines) == len(scored) + 1


def test_snp_file_parser_round_trips_the_extract_output(tmp_path):
    path = str(tmp_path / "snps.txt")
    with open(path, "w") as handle:
        handle.write("# a comment\n\n100 A.1 G C\n200 None T A\n")
    assert snp_file_parser(path) == [[100, "A.1", "G", "C"], [200, "None", "T", "A"]]
