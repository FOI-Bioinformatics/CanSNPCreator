"""Behaviour that is expected to survive the hardening work.

These tests describe what the tools do today and are intended to keep doing:
the parsers, the five score terms and their scaling, and the output columns.
Tests for the defects being fixed live beside them, one file per area.
"""

import pytest
from conftest import GENOMES, INCONGRUENT, PLANTED, node_isolates

from cansnpscorer.helpers.data_utils import (
    apply_filters,
    codon_comparer,
    codon_complement,
    gc_content,
    neighbour_distance,
    search,
    snp_extractor,
    summarize_scores,
)
from cansnpscorer.helpers.file_utils import read_alignment, read_gff3, read_intervals
from cansnpscorer.helpers.genetic_code import codon_table
from cansnpscorer.helpers.names import clean_name
from cansnpscorer.helpers.scoring_methods import MAXIMUM_SCORE, TERMS, ScoredSnp

# -- parsers -----------------------------------------------------------------


def test_the_alignment_is_read_in_file_order(dataset):
    alignment = read_alignment(dataset.alignment)
    assert list(alignment) == GENOMES
    assert alignment["ref"] == dataset.sequences["ref"]


def test_clean_name_writes_every_separator_alike():
    assert clean_name("GCA_000007125.1 ASM712v1") == "GCA_000007125_1_ASM712v1"
    assert clean_name("S1_2") != clean_name("S12")


def test_gff3_genes_and_cds_are_read_per_sequence_zero_based(dataset):
    model = read_gff3(dataset.gff)
    assert sorted(model.genes) == ["c1", "c2"]
    assert model.genes["c1"] == [(150, 449), (550, 649)], "gene and CDS once"
    assert model.genes["c2"] == [(10, 199)]
    assert model.cds["c1"] == [(150, 449, "+", 0), (550, 649, "-", 0)]


def test_bed_regions_are_read_per_sequence_zero_based(dataset):
    # BED is 0-based and half open; the regions are held 0-based and inclusive,
    # as the GFF3 genes are.
    assert read_intervals(dataset.repeats) == {"c1": [(59, 94)], "c2": [(119, 159)]}


# -- extraction --------------------------------------------------------------


def test_snp_extractor_finds_the_biallelic_columns(dataset):
    snps = snp_extractor(read_alignment(dataset.alignment))
    found = [snp["index"] for snp in snps]
    expected = sorted([column for column, _, _ in PLANTED.values()] + [INCONGRUENT[0]])
    assert found == expected, "the tri-allelic and gapped columns must be left out"


def test_cansnp_extractor_assigns_every_planted_node(cansnps):
    assigned = {cansnp.column: cansnp.node for cansnp in cansnps}
    for node, (column, _ancestral, _derived) in PLANTED.items():
        assert assigned[column] == node
    assert assigned[INCONGRUENT[0]] is None, "a SNP crossing two branches is not canonical"


def test_cansnp_extractor_records_the_polarity(cansnps):
    # The base shared below the node is the derived one, which the extractor
    # used to throw away by reporting the two bases in an arbitrary order.
    planted = {column: (ancestral, derived) for column, ancestral, derived in PLANTED.values()}
    for cansnp in cansnps:
        if cansnp.column not in planted:
            continue
        ancestral, derived = planted[cansnp.column]
        assert (cansnp.derived, cansnp.ancestral) == (derived, ancestral), cansnp.node


def test_hierarchy_gives_internal_nodes_their_isolates(isolates_by_node):
    for node in ["root", "A", "B", "A.1", "A.2", "B.1", "B.2"]:
        assert sorted(isolates_by_node[node]) == sorted(node_isolates(node)), node


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


def test_failing_scores_are_marked_and_kept():
    scores = [make_score(1, "A", 300), make_score(2, "A", 100), make_score(3, "B", 400)]
    apply_filters(scores, 150)
    assert [entry.filters for entry in scores] == [[], ["LowScore"], []]

    scores = [make_score(1, "A", 300), make_score(2, "A", 200), make_score(3, "B", 400)]
    apply_filters(scores, 50, top=1)
    assert [entry.filters for entry in scores] == [[], ["NotTop"], []]


def test_the_best_of_a_node_are_taken_apart_from_each_other():
    # 300 and 290 lie 10 bases apart, so with a spacing of 100 the second best
    # gives way to the third, which lies far from both.
    scores = [make_score(1000, "A", 300), make_score(1010, "A", 290),
              make_score(5000, "A", 280), make_score(1005, "B", 250)]
    apply_filters(scores, 50, top=2, spacing=100)
    assert [entry.filters for entry in scores] == [[], ["NotTop"], [], []]
    # Without spacing the two best are taken, however close.
    for entry in scores:
        entry.filters.clear()
    apply_filters(scores, 50, top=2)
    assert [entry.filters for entry in scores] == [[], [], ["NotTop"], []]


def test_ties_go_to_the_earlier_position_as_export_breaks_them():
    scores = [make_score(3, "A", 300), make_score(1, "A", 300)]
    apply_filters(scores, 50, top=1)
    assert [entry.filters for entry in scores] == [["NotTop"], []]


def test_summarize_scores_counts_snps_per_node():
    scores = [make_score(1, "A", 300), make_score(2, "A", 300), make_score(3, "B", 300)]
    assert summarize_scores(["A", "B", "C"], scores) == {"A": 2, "B": 1}
