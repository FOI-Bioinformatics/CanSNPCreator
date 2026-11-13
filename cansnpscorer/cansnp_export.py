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

A scheme may also carry more than one canSNP per node. CanSNPer2 calls one SNP
per node, so the SNP file keeps the best of each, and the full set goes to an
optional <organism>_canSNPs.txt, ranked within each node, whose rank 1 is the
SNP of the SNP file.
"""

import logging
import os
import re
from collections.abc import Sequence

from cansnpscorer import vcf
from cansnpscorer.helpers.data_utils import ranking, select_per_node
from cansnpscorer.helpers.file_utils import (
    fasta_ids,
    read_genome_map,
    read_hierarchy,
    reference_parser,
)
from cansnpscorer.helpers.scoring_methods import ScoredSnp
from cansnpscorer.hierarchy import Hierarchy, HierarchyError
from cansnpscorer.model import Reference, ReferenceError
from cansnpscorer.validate import report

GENBANK = re.compile(r"GCA_\d{9}\.\d+")
REFSEQ = re.compile(r"GCF_\d{9}\.\d+")
#: An organism names the files, so it is limited to what a file name holds safely.
ORGANISM = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")
BASES = frozenset("ACGT")

TREE_SUFFIX = "_tree.txt"
SNP_SUFFIX = "_snp.txt"
REFERENCES_SUFFIX = "_references.txt"
GENOMES_SUFFIX = "_genomes.txt"
SNP_SET_SUFFIX = "_canSNPs.txt"
SNP_SET_COLUMNS = ("node", "rank", "genome", "position", "derived_base", "ancestral_base",
                   "score")
GENOME_COLUMNS = ("genome", "node", "role")
ACCESSION = re.compile(r"GC[AF]_\d{9}\.\d+")

SNP_COLUMNS = ("snp_id", "strain", "reference", "genome", "position",
               "derived_base", "ancestral_base")
REFERENCE_COLUMNS = ("genome", "strain", "genbank_id", "refseq_id", "assembly_name")


def best_per_node(scores: Sequence[ScoredSnp]) -> dict[str, ScoredSnp]:
    """The highest scoring SNP of each node, among those given.

    A scheme carries exactly one SNP per node, and the caller passes only the
    canSNPs that passed evaluate's filters. Ties are broken by contig and
    position, so that the same input always exports the same scheme.
    """
    best: dict[str, ScoredSnp] = {}
    for score in scores:
        if score.node not in best or ranking(score) < ranking(best[score.node]):
            best[score.node] = score
    return best


def check_metadata(genome: str, strain: str, reference: str, genbank_id: str,
                   refseq_id: str, assembly_name: str, organism: str = "x") -> list[str]:
    """The reference genome's row must satisfy the CanSNPer2-data rules."""
    problems: list[str] = []
    if not ORGANISM.fullmatch(organism or ""):
        problems.append(
            f"the organism '{organism}' names the files, so it may hold only letters, digits, "
            f"'_', '.' and '-', as francisella_tularensis does")

    if "_" in genome:
        problems.append(
            f"the genome name '{genome}' holds an underscore; CanSNPer2 recovers the reference "
            f"name from an alignment file name by splitting on the first underscore")
    if not genome:
        problems.append("the genome name is empty")
    if not GENBANK.fullmatch(genbank_id):
        problems.append(
            f"'{genbank_id}' is not a well formed GenBank accession, such as GCA_000007125.1")
    if not REFSEQ.fullmatch(refseq_id):
        problems.append(
            f"'{refseq_id}' is not a well formed RefSeq accession, such as GCF_000007125.1")
    if not assembly_name:
        problems.append("the assembly name is empty")

    for label, value in [("genome", genome), ("strain", strain), ("reference", reference),
                         ("assembly name", assembly_name)]:
        if ";" in value:
            problems.append(
                f"the {label} '{value}' holds a semicolon, which the tree reader cannot carry")
        if any(character in value for character in "\t\r\n"):
            problems.append(
                f"the {label} {value!r} holds a tab or a line break, which would shift the "
                f"columns of the files")
        if value != value.strip():
            problems.append(f"the {label} '{value}' has leading or trailing whitespace")

    return problems


def check_scheme(hierarchy: Hierarchy, best: dict[str, ScoredSnp],
                 reference: Reference) -> tuple[list[str], list[str]]:
    """Check the scheme against the rules of the CanSNPer2-data source files."""
    problems: list[str] = []
    warnings: list[str] = []

    for node in hierarchy.nodes:
        if ";" in node:
            problems.append(f"node '{node}' holds a semicolon, which the tree reader cannot carry")

    unknown = sorted(set(best) - set(hierarchy.parents))
    if unknown:
        problems.append(
            f"the scores name {len(unknown)} nodes that the hierarchy does not hold: "
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


def paths_for(hierarchy: Hierarchy, kept: Sequence[str]) -> list[list[str]]:
    """One path per kept node, from the root, skipping the nodes left out."""
    keeping = set(kept)
    paths = []
    for node in hierarchy.nodes:
        if node not in keeping:
            continue
        paths.append([name for name in hierarchy.path(node) if name in keeping])
    return paths


def genome_rows(hierarchy: Hierarchy, kept: Sequence[str], clusters: dict[str, list[str]],
                aligned: Sequence[str]) -> tuple[list[tuple[str, str, str]], list[str]]:
    """The rows of the genomes file, and warnings about genomes left out.

    A genome placed at a node the scheme left out moves to the nearest node
    above it that was kept, as that node's children do. A genome whose
    sequence is in the alignment is a representative; one placed only through
    the genome map is a member; an aligned genome placed at no node is the
    outgroup. Each is written as the assembly accession its name holds.
    """
    from cansnpscorer.helpers.names import resolve_names

    keeping = set(kept)

    def kept_at(node: str) -> str | None:
        while node is not None and node not in keeping:
            node = hierarchy.parents.get(node)
        return node

    placed = {genome: node for node, genomes in clusters.items() for genome in genomes}
    matched, _unmatched, _ambiguous = resolve_names(list(placed), aligned)
    aligned_placed = set(matched.values())

    rows: list[tuple[str, str, str]] = []
    warnings: list[str] = []
    no_accession: list[str] = []
    seen: dict[str, str] = {}

    def add(name: str, node: str, role: str) -> None:
        found = ACCESSION.search(name)
        if found is None:
            no_accession.append(name)
            return
        accession = found.group(0)
        assembly = accession[4:].split(".")[0]
        if assembly in seen:
            return
        seen[assembly] = accession
        rows.append((accession, node, role))

    for genome, node in placed.items():
        if node not in hierarchy:
            continue
        add(genome, kept_at(node) or hierarchy.root,
            "representative" if genome in matched else "member")
    for record in aligned:
        if record not in aligned_placed:
            add(record, "", "outgroup")

    if no_accession:
        warnings.append(
            f"{len(no_accession)} genomes hold no assembly accession in their names and were "
            f"left out of the genomes file: {', '.join(sorted(no_accession)[:5])}")
    order = {"representative": 0, "member": 1, "outgroup": 2}
    rows.sort(key=lambda row: (order[row[2]], row[1], row[0]))
    return rows, warnings


def write_genomes(outdir: str, organism: str, rows: Sequence[tuple[str, str, str]]) -> str:
    path = os.path.join(outdir, organism + GENOMES_SUFFIX)
    with open(path, "w", newline="\n") as handle:
        handle.write("\t".join(GENOME_COLUMNS) + "\n")
        handle.writelines("\t".join(row) + "\n" for row in rows)
    return path


def write_snp_set(outdir: str, organism: str, paths: Sequence[Sequence[str]],
                  chosen: dict[str, Sequence[ScoredSnp]], reference: Reference,
                  genome: str) -> str:
    """Write every chosen canSNP of each kept node, ranked within the node.

    Positions are over the concatenated genome, as in the SNP file, and rank 1
    of each node is the SNP the SNP file carries for it.
    """
    path = os.path.join(outdir, organism + SNP_SET_SUFFIX)
    with open(path, "w", newline="\n") as handle:
        handle.write("\t".join(SNP_SET_COLUMNS) + "\n")
        for tree_path in paths:
            node = tree_path[-1]
            for rank, score in enumerate(chosen[node], start=1):
                position = reference.contig(score.contig).start + score.position
                handle.write("\t".join([
                    node, str(rank), genome, str(position), score.derived, score.ancestral,
                    str(score.total),
                ]) + "\n")
    return path


def write_scheme(outdir: str, organism: str, paths: Sequence[Sequence[str]],
                 best: dict[str, ScoredSnp], reference: Reference, genome: str, strain: str,
                 citation: str, genbank_id: str, refseq_id: str,
                 assembly_name: str) -> list[str]:
    """Write the three source files and return their paths."""
    os.makedirs(outdir, exist_ok=True)
    written = []

    tree_path = os.path.join(outdir, organism + TREE_SUFFIX)
    with open(tree_path, "w", newline="\n") as handle:
        handle.writelines("\t".join(path) + "\n" for path in paths)
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
         assembly_name=None, prune=False, genome_map=None, alignment=None,
         snps_per_node=1, spacing=1000) -> int:
    """Write a scored scheme as CanSNPer2-data source files. Returns an exit status.

    ``snps_per_node`` is how many canSNPs of each node to keep, at least
    ``spacing`` bases apart, or None for all of them. With more than one the
    full set is written to <organism>_canSNPs.txt.
    """
    if bool(genome_map) != bool(alignment):
        logging.error("the genomes file needs both --genome-map and --alignment: where each "
                      "genome sits, and which genomes were aligned")
        return 1
    try:
        clusters = read_genome_map(genome_map) if genome_map else None
        aligned = fasta_ids(alignment) if alignment else None
        tree = read_hierarchy(hierarchy)
        sequence = reference_parser(reference)
        scored = [score for score in vcf.read_scores(scores, sequence) if score.passed]
    except (OSError, ValueError, HierarchyError, ReferenceError) as error:
        logging.error(error)
        return 1

    chosen = select_per_node(scored, snps_per_node, spacing)
    best = {node: picked[0] for node, picked in chosen.items() if picked}
    problems, warnings = check_scheme(tree, best, sequence)
    problems += check_metadata(genome, strain, citation, genbank_id, refseq_id,
                               assembly_name, organism)

    without = [node for node in tree.nodes if node not in best]
    if without:
        if tree.root in without:
            problems.append(
                f"the root '{tree.root}' has no canSNP, and a scheme cannot be written "
                f"without one. A root's canSNP needs a genome outside the root in the "
                f"alignment, an outgroup placed at no node (collapse --outgroup); with one, "
                f"lower the evaluate cutoff, or root the scheme elsewhere")
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
    if snps_per_node != 1:
        written.append(write_snp_set(outdir, organism, paths, chosen, sequence, genome))
        counts = sorted(len(chosen[path[-1]]) for path in paths)
        short = [path[-1] for path in paths
                 if snps_per_node is not None and len(chosen[path[-1]]) < snps_per_node]
        logging.info(
            f"The canSNP set holds {sum(counts)} canSNPs, {counts[0]} to {counts[-1]} per node")
        if short:
            logging.warning(
                f"{len(short)} nodes have fewer than {snps_per_node} canSNPs that passed and lie "
                f"{spacing} bases apart: {', '.join(short[:5])}")
    if clusters is not None:
        rows, notes = genome_rows(tree, kept, clusters, aligned)
        for note in notes:
            logging.warning(note)
        written.append(write_genomes(outdir, organism, rows))
        counts = {role: sum(1 for row in rows if row[2] == role)
                  for role in ("representative", "member", "outgroup")}
        logging.info("The genomes file lists " + ", ".join(
            f"{count} {role}s" for role, count in counts.items()))

    logging.info(f"Wrote a scheme of {len(kept)} nodes: {', '.join(written)}")
    logging.info(
        "Check it with the CanSNPer2-data repository's validate.py before building a database")
    return 0
