"""Regression tests for the third audit round: rooting on outgroups, nested
nodes under --min-called, alleles checked on export, and stale export files."""

import logging
import os

import pytest
from conftest import CLUSTERS, vcf_rewrite
from test_export import METADATA, cli_run, several  # noqa: F401  (several is a fixture)
from test_export import scheme as scheme

from cansnpscorer import vcf
from cansnpscorer.cansnp_export import check_scheme, genome_rows
from cansnpscorer.cansnp_export import main as export_main
from cansnpscorer.helpers.data_utils import canSNP_extractor, snp_extractor
from cansnpscorer.helpers.scoring_methods import TERMS, ScoredSnp, contexts_for
from cansnpscorer.hierarchy import Hierarchy, HierarchyError
from cansnpscorer.model import Reference
from cansnpscorer.trees.collapsing import collapse
from cansnpscorer.trees.rooting import reroot


def clades(tree: Hierarchy) -> set[frozenset[str]]:
    return {frozenset(tree.leaves_below(node)) for node in tree.nodes if tree.children[node]}


# -- rooting on the outgroups ------------------------------------------------


def test_a_tree_rooted_inside_the_ingroup_is_rooted_on_the_outgroup():
    tree = Hierarchy.from_newick("(((O:1,A:1)90:2,B:1)80:3,(C:1,D:1)70:4);", support_labels=True)
    rooted = reroot(tree, ["O"])
    assert rooted.parents["O"] == rooted.root
    assert clades(rooted) == {frozenset("OABCD"), frozenset("ABCD"), frozenset("BCD"), frozenset("CD")}
    # Each support stays with its split of the genomes.
    bcd = next(n for n in rooted.nodes if set(rooted.leaves_below(n)) == set("BCD"))
    cd = next(n for n in rooted.nodes if set(rooted.leaves_below(n)) == set("CD"))
    assert rooted.support[bcd] == pytest.approx(0.9)
    assert rooted.support[cd] == pytest.approx(0.7)
    # The old root had two children; its two branches are now one.
    assert rooted.lengths[cd] == pytest.approx(7)
    # The branch the root was placed on is halved.
    assert rooted.lengths["O"] == pytest.approx(0.5)


def test_a_tree_already_rooted_on_the_outgroup_is_unchanged():
    for text in ("(O,(A,(B,(C,D))));", "(A,O,(B,C));", "((O,P),(A,B));"):
        tree = Hierarchy.from_newick(text, support_labels=True)
        outgroups = ["O", "P"] if "P" in text else ["O"]
        assert reroot(tree, outgroups) is tree


def test_several_outgroups_are_rooted_on_together():
    tree = Hierarchy.from_newick("((A,B),((O,P),C));", support_labels=True)
    rooted = reroot(tree, ["O", "P"])
    assert frozenset("OP") in clades(rooted)
    assert frozenset("ABC") in clades(rooted)


def test_outgroups_that_form_no_clade_are_refused():
    with pytest.raises(HierarchyError, match="do not form a clade"):
        reroot(Hierarchy.from_newick("((O,A),(P,B));"), ["O", "P"])


def test_collapse_finds_the_clades_of_the_tree_rooted_on_the_outgroup():
    result = collapse(Hierarchy.from_newick("(((O,A),B),(C,D));", support_labels=True), outgroups=["O"])
    placed = {}
    for genome, node in result.genomes.items():
        placed.setdefault(node, set()).add(genome)
    below = {node: set() for node in result.hierarchy.nodes}
    for node in result.hierarchy.nodes:
        for genome, at in result.genomes.items():
            if at in [node, *result.hierarchy.descendants(node)]:
                below[node].add(genome)
    assert sorted(map(sorted, below.values())) == [["A", "B", "C", "D"], ["B", "C", "D"], ["C", "D"]]
    assert "O" not in result.genomes


# -- nested nodes under --min-called ------------------------------------------


def test_a_snp_that_fits_a_node_and_its_child_is_not_assigned(caplog):
    # R -> P -> C, and R -> Q. c is the only genome telling P from C, and has
    # no base at the SNP.
    alignment = {"ref": "CA", "O": "CA", "a": "TA", "b": "TA", "c": "-A", "d": "CC", "e": "CC"}
    groups = {"R": ["a", "b", "c", "d", "e"], "P": ["a", "b", "c"], "C": ["a", "b"], "Q": ["d", "e"]}
    snps = snp_extractor(alignment, min_called=0.6)
    with caplog.at_level(logging.WARNING):
        found = {snp.column: snp for snp in canSNP_extractor(alignment, groups, snps, 0.6)}
    assert found[0].node is None
    assert "fit both a node and a node below it" in caplog.text
    assert "outgroup genome" not in caplog.text
    # With c called, the column marks C alone.
    alignment["c"] = "CA"
    snps = snp_extractor(alignment, min_called=0.6)
    found = {snp.column: snp for snp in canSNP_extractor(alignment, groups, snps, 0.6)}
    assert found[0].node == "C"


# -- alleles checked on export -------------------------------------------------


def test_scores_whose_ancestral_allele_is_neither_ref_nor_alt_are_refused(dataset, scheme, tmp_path):
    scores, outdir = scheme

    def third_base(fields):
        if "AA=" not in fields[7]:
            return fields
        third = next(base for base in "ACGT" if base not in (fields[3], fields[4]))
        fields[7] = ";".join(
            f"AA={third}" if item.startswith("AA=") else item for item in fields[7].split(";")
        )
        return fields

    wrong = vcf_rewrite(scores, tmp_path / "wrong.vcf", change=third_base)
    with pytest.raises(vcf.VcfError, match="neither REF nor ALT"):
        vcf.read_scores(wrong)
    assert export_main(wrong, dataset.hierarchy, dataset.reference, outdir, **METADATA) == 1
    assert not os.path.exists(outdir)


def test_a_snp_whose_alleles_miss_the_reference_base_is_a_problem():
    reference = Reference.from_records([("c1", "ACGTACGTAC")])
    tree = Hierarchy([("root", None)])
    best = {
        "root": ScoredSnp(
            contig="c1",
            position=1,
            node="root",
            total=500,
            terms=dict.fromkeys(TERMS, 100),
            derived="G",
            ancestral="T",
        )
    }
    problems, _warnings = check_scheme(tree, best, reference)
    assert any("reference base there is A" in problem for problem in problems)
    best["root"].ancestral = "A"
    assert check_scheme(tree, best, reference) == ([], [])


# -- the genomes file and files of an earlier export --------------------------


def test_genomes_at_a_node_the_hierarchy_lacks_are_reported():
    tree = Hierarchy([("root", None), ("A", "root")])
    rows, warnings = genome_rows(
        tree,
        ["root", "A"],
        {"A": ["GCF_000000001.1"], "Z": ["GCF_000000002.1"]},
        ["GCF_000000001.1", "GCF_000000002.1"],
    )
    assert [row[0] for row in rows] == ["GCF_000000001.1"]
    assert any("GCF_000000002.1 (Z)" in warning for warning in warnings)


def test_files_an_earlier_export_wrote_are_removed(dataset, several):  # noqa: F811
    scores, outdir = several
    organism = METADATA["organism"]
    assert (
        export_main(
            scores,
            dataset.hierarchy,
            dataset.reference,
            outdir,
            snps_per_node=3,
            genome_map=dataset.genome_map,
            alignment=dataset.alignment,
            **METADATA,
        )
        == 0
    )
    cansnp_file = os.path.join(outdir, organism + "_canSNPs.txt")
    old_snp_file = os.path.join(outdir, organism + "_snp.txt")
    genomes_file = os.path.join(outdir, organism + "_genomes.txt")
    assert os.path.exists(cansnp_file) and os.path.exists(genomes_file)
    # The SNP file of an export before the canSNP file.
    with open(old_snp_file, "w") as handle:
        handle.write("snp_id\tstrain\n")
    assert export_main(scores, dataset.hierarchy, dataset.reference, outdir, **METADATA) == 0
    assert not os.path.exists(old_snp_file) and not os.path.exists(genomes_file)
    assert os.path.exists(cansnp_file), "the canSNP file is written at every export"


# -- evaluate and extract ------------------------------------------------------


def test_contexts_share_the_interval_lists_of_their_contig():
    reference = Reference.from_records([("c1", "ACGT" * 50)])
    contig = reference.contig("c1")
    genes = {"c1": [(100, 150), (10, 60)]}
    contexts = contexts_for([(contig, 20, "A", ("A", "C")), (contig, 120, "A", ("A", "C"))], genes, {}, {})
    assert contexts[0].genes is contexts[1].genes
    assert contexts[0].genes == [(10, 60), (100, 150)]


def test_extract_reports_the_number_of_genomes_of_the_scheme(dataset, tmp_path, caplog):
    with caplog.at_level(logging.INFO):
        assert (
            cli_run(
                "extract",
                "--genome-map",
                dataset.genome_map,
                "--hierarchy",
                dataset.hierarchy,
                "--alignment",
                dataset.alignment,
                "--reference",
                dataset.reference,
                "-o",
                str(tmp_path / "s.vcf"),
            )
            == 0
        )
    genomes = len({name for names in CLUSTERS.values() for name in names})
    assert f"among the {genomes} genomes of the scheme" in caplog.text
