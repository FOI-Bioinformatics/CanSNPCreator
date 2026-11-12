"""Export of a scored scheme as CanSNPer2-data source files.

CanSNPer2-data holds a scheme as three tab separated files: the tree as one
line per path from the root, one SNP per node, and the reference genomes. The
rules those files must follow are set out in that repository's source/README.md
and enforced by its validate.py, and the checks here are the same ones, so that
what this writes can be added to that repository without hand editing:

- every node of the tree carries exactly one SNP
- positions are 1-based over the concatenated reference genome, since the SNP
  file has no contig column
- a genome name holds no underscore, because CanSNPer2 recovers the reference
  name from an alignment file name by splitting on the first one
- no semicolon appears anywhere, because the tree reader turns tabs into
  semicolons before splitting
- accessions are well formed, because CanSNPer2-download builds NCBI paths out
  of their digits
"""

import logging
import os
import re
from typing import Dict, List, Sequence, Tuple

from cansnpscorer import vcf
from cansnpscorer.helpers.file_utils import read_hierarchy, reference_parser
from cansnpscorer.helpers.scoring_methods import ScoredSnp
from cansnpscorer.hierarchy import Hierarchy, HierarchyError
from cansnpscorer.model import Reference, ReferenceError
from cansnpscorer.validate import report

GENBANK = re.compile(r"^GCA_\d{9}\.\d+$")
REFSEQ = re.compile(r"^GCF_\d{9}\.\d+$")
BASES = frozenset("ACGT")

TREE_SUFFIX = "_tree.txt"
SNP_SUFFIX = "_snp.txt"
REFERENCES_SUFFIX = "_references.txt"

SNP_COLUMNS = ("snp_id", "strain", "reference", "genome", "position",
               "derived_base", "ancestral_base")
REFERENCE_COLUMNS = ("genome", "strain", "genbank_id", "refseq_id", "assembly_name")


def best_per_node(scores: Sequence[ScoredSnp]) -> Dict[str, ScoredSnp]:
    """The highest scoring SNP of each node, among those given.

    A scheme carries exactly one SNP per node, and the caller passes only the
    canSNPs that passed evaluate's filters. Ties are broken by contig and
    position, so that the same input always exports the same scheme.
    """
    def ranking(score: ScoredSnp):
        return (-score.total, score.contig, score.position)

    best: Dict[str, ScoredSnp] = {}
    for score in scores:
        if score.node not in best or ranking(score) < ranking(best[score.node]):
            best[score.node] = score
    return best


def check_metadata(genome: str, strain: str, reference: str, genbank_id: str,
                   refseq_id: str, assembly_name: str) -> List[str]:
    """The reference genome's row must satisfy the CanSNPer2-data rules."""
    problems: List[str] = []

    if "_" in genome:
        problems.append(
            f"the genome name '{genome}' holds an underscore; CanSNPer2 recovers the reference "
            f"name from an alignment file name by splitting on the first underscore")
    if not genome:
        problems.append("the genome name is empty")
    if not GENBANK.match(genbank_id):
        problems.append(
            f"'{genbank_id}' is not a well formed GenBank accession, such as GCA_000007125.1")
    if not REFSEQ.match(refseq_id):
        problems.append(
            f"'{refseq_id}' is not a well formed RefSeq accession, such as GCF_000007125.1")
    if not assembly_name:
        problems.append("the assembly name is empty")

    for label, value in [("genome", genome), ("strain", strain), ("reference", reference),
                         ("assembly name", assembly_name)]:
        if ";" in value:
            problems.append(
                f"the {label} '{value}' holds a semicolon, which the tree reader cannot carry")
        if value != value.strip():
            problems.append(f"the {label} '{value}' has leading or trailing whitespace")

    return problems


def check_scheme(hierarchy: Hierarchy, best: Dict[str, ScoredSnp],
                 reference: Reference) -> Tuple[List[str], List[str]]:
    """Check the scheme against the rules of the CanSNPer2-data source files."""
    problems: List[str] = []
    warnings: List[str] = []

    for node in hierarchy.nodes:
        if ";" in node:
            problems.append(f"node '{node}' holds a semicolon, which the tree reader cannot carry")

    unknown = sorted(set(best) - set(hierarchy.parents))
    if unknown:
        problems.append(
            f"the score file names {len(unknown)} nodes that the hierarchy does not hold: "
            f"{', '.join(unknown[:5])}")

    for node, score in sorted(best.items()):
        if node not in hierarchy:
            continue
        if score.contig not in reference.names:
            problems.append(
                f"the SNP for '{node}' is on contig '{score.contig}', which the reference does "
                f"not hold")
            continue
        contig = reference.contig(score.contig)
        if not 1 <= score.position <= contig.length:
            problems.append(
                f"the SNP for '{node}' is at {score.contig}:{score.position}, outside that "
                f"contig of {contig.length} bases")
        if not score.derived or not score.ancestral:
            problems.append(f"the SNP for '{node}' records no derived and ancestral base")
        elif score.derived == score.ancestral:
            problems.append(
                f"the SNP for '{node}' gives the same base, {score.derived}, as derived and "
                f"ancestral")
        elif not {score.derived, score.ancestral} <= BASES:
            problems.append(
                f"the SNP for '{node}' uses bases outside A, C, G and T: "
                f"{score.derived} and {score.ancestral}")

    return problems, warnings


def paths_for(hierarchy: Hierarchy, kept: Sequence[str]) -> List[List[str]]:
    """One path per kept node, from the root, skipping the nodes left out."""
    keeping = set(kept)
    paths = []
    for node in hierarchy.nodes:
        if node not in keeping:
            continue
        paths.append([name for name in hierarchy.path(node) if name in keeping])
    return paths


def write_scheme(outdir: str, organism: str, paths: Sequence[Sequence[str]],
                 best: Dict[str, ScoredSnp], reference: Reference, genome: str, strain: str,
                 citation: str, genbank_id: str, refseq_id: str,
                 assembly_name: str) -> List[str]:
    """Write the three source files and return their paths."""
    os.makedirs(outdir, exist_ok=True)
    written = []

    tree_path = os.path.join(outdir, organism + TREE_SUFFIX)
    with open(tree_path, "w", newline="\n") as handle:
        for path in paths:
            handle.write("\t".join(path) + "\n")
    written.append(tree_path)

    snp_path = os.path.join(outdir, organism + SNP_SUFFIX)
    with open(snp_path, "w", newline="\n") as handle:
        handle.write("\t".join(SNP_COLUMNS) + "\n")
        for path in paths:
            node = path[-1]
            score = best[node]
            # The SNP file has no contig column, so positions are counted over
            # the concatenated genome, as the published multi-contig schemes do.
            position = reference.contig(score.contig).start + score.position
            handle.write("\t".join([
                node, strain, citation, genome, str(position),
                score.derived, score.ancestral,
            ]) + "\n")
    written.append(snp_path)

    references_path = os.path.join(outdir, organism + REFERENCES_SUFFIX)
    with open(references_path, "w", newline="\n") as handle:
        handle.write("\t".join(REFERENCE_COLUMNS) + "\n")
        handle.write("\t".join([genome, strain, genbank_id, refseq_id, assembly_name]) + "\n")
    written.append(references_path)

    return written


def main(scores=None, hierarchy=None, reference=None, outdir=".", organism=None,
         genome=None, strain=None, citation=None, genbank_id=None, refseq_id=None,
         assembly_name=None, prune=False) -> int:
    """Write a scored scheme as CanSNPer2-data source files. Returns an exit status."""
    try:
        tree = read_hierarchy(hierarchy)
        sequence = reference_parser(reference)
        scored = [score for score in vcf.read_scores(scores, sequence) if score.passed]
    except (OSError, ValueError, HierarchyError, ReferenceError) as error:
        logging.error(error)
        return 1

    best = best_per_node(scored)
    problems, warnings = check_scheme(tree, best, sequence)
    problems += check_metadata(genome, strain, citation, genbank_id, refseq_id,
                               assembly_name)

    without = [node for node in tree.nodes if node not in best]
    if without:
        if tree.root in without:
            problems.append(
                f"the root '{tree.root}' has no canSNP, and a scheme cannot be written "
                f"without one. Lower the evaluate cutoff, or root the scheme elsewhere")
        if prune:
            warnings.append(
                f"{len(without)} nodes have no canSNP and were left out of the scheme, so their "
                f"children now hang from the nearest node that has one: {', '.join(without[:5])}")
        else:
            problems.append(
                f"{len(without)} nodes have no canSNP: {', '.join(without[:5])}. Every node of a "
                f"CanSNPer2 scheme carries exactly one, so either lower the evaluate cutoff or pass "
                f"--prune to leave these nodes out")

    if not report(problems, warnings, "scheme"):
        return 2

    kept = [node for node in tree.nodes if node in best]
    paths = paths_for(tree, kept)
    written = write_scheme(outdir, organism, paths, best, sequence, genome, strain,
                           citation, genbank_id, refseq_id, assembly_name)

    logging.info(f"Wrote a scheme of {len(kept)} nodes: {', '.join(written)}")
    logging.info(
        "Check it with the CanSNPer2-data repository's validate.py before building a database")
    return 0
