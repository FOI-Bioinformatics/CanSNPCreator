"""Export of a scored scheme as CanSNPer2-data source files.

CanSNPer2-data holds a scheme as three tab separated files: the tree as one
line per path from the root, one SNP per node, and the reference genomes. The
rules those files must follow are set out in that repository's schemes/README.md
and enforced by its validate.py, and the checks here are the same ones, so that
what this writes can be added to that repository without hand editing:

- every node of the tree carries exactly one SNP
- each SNP names its reference sequence and gives a 1-based position on that
  sequence, as CanSNPer2 2.3.0 reads them
- a sequence is named by its GenBank accession, as in the FASTA files
  CanSNPer2-download fetches. A sequence the reference names by its RefSeq
  accession (NC_, NZ_ or AC_) is written under the GenBank accession that
  --contig-aliases pairs it with, and the export stops when there is none. A
  name that is neither is written as it is, with a warning
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
    read_contig_aliases,
    read_genome_map,
    read_hierarchy,
    reference_parser,
)
from cansnpscorer.helpers.scoring_methods import ScoredSnp
from cansnpscorer.hierarchy import Hierarchy, HierarchyError
from cansnpscorer.model import Reference, ReferenceError
from cansnpscorer.scheme import CANSNP_COLUMNS, CANSNP_SUFFIX
from cansnpscorer.validate import report

#: The version of an accession has one or two digits, as CanSNPer2-data's
#: validate.py requires; a longer one is refused here rather than there. The
#: names module matches any version, since it only finds a genome by its name
#: and writes nothing.
VERSION = r"\.\d{1,2}(?!\d)"
GENBANK = re.compile(r"GCA_\d{9}" + VERSION)
REFSEQ = re.compile(r"GCF_\d{9}" + VERSION)
#: An organism names the files, so it is limited to what a file name holds safely.
ORGANISM = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")
BASES = frozenset("ACGT")

TREE_SUFFIX = "_tree.txt"
SNP_SUFFIX = "_snp.txt"
REFERENCES_SUFFIX = "_references.txt"
GENOMES_SUFFIX = "_genomes.txt"

GENOME_COLUMNS = ("genome", "node", "role")
ACCESSION = re.compile(r"GC[AF]_\d{9}" + VERSION)
#: An accession in any version or none, to say why a name gave no accession.
ANY_ACCESSION = re.compile(r"GC[AF]_\d{9}")
#: A sequence named by its RefSeq accession, as NC_002971.4.
REFSEQ_SEQUENCE = re.compile(r"(NC|NZ|AC)_")
#: A RefSeq contig or scaffold, as NT_ or NW_, which is written as it is.
OTHER_REFSEQ = re.compile(r"(NT|NW)_")
#: A sequence named by its GenBank (INSDC) accession with a version, as
#: AE016828.3 or JAABCD010000001.1: letters, then digits, then the version.
GENBANK_SEQUENCE = re.compile(r"[A-Z]{1,6}[0-9]+\.[0-9]+")

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
            f"'{genbank_id}' is not a well formed GenBank accession, such as GCA_000007125.1, "
            f"with a version of one or two digits")
    if not REFSEQ.fullmatch(refseq_id):
        problems.append(
            f"'{refseq_id}' is not a well formed RefSeq accession, such as GCF_000007125.1, "
            f"with a version of one or two digits")
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
        # Every file of the scheme is tab separated with one row per line, so
        # such a name would not read back. A quoted Newick label can hold one.
        if any(character in node for character in "\t\r\n"):
            problems.append(
                f"node {node!r} holds a tab or a line break, which would shift the columns "
                f"and lines of the files")

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
        elif 1 <= score.position <= contig.length:
            base = contig.base(score.position - 1)
            if base in BASES and base not in (score.derived, score.ancestral):
                problems.append(
                    f"the SNP for '{node}' at {score.contig}:{score.position} gives "
                    f"{score.derived}/{score.ancestral}, but the reference base there is {base}")

    return problems, warnings


def sequence_names(contigs: Sequence[str], aliases: Sequence[set[str]] | None,
                   alias_file: str | None = None,
                   ) -> tuple[dict[str, str], list[str], list[str]]:
    """The name each reference sequence is written under in the scheme.

    CanSNPer2 2.3.0 finds a canSNP's sequence by name in the reference FASTA
    that CanSNPer2-download fetches from GenBank, so a sequence is written
    under its GenBank accession. The rule, for each sequence name:

    - a GenBank accession with a version is written as it is;
    - a RefSeq accession (NC_, NZ_ or AC_) is written as the one name of its
      alias set, from --contig-aliases, that is a GenBank accession with a
      version and no underscore. With no alias table, no such name, or more
      than one, the export stops. It also stops when the alias table, an
      NCBI assembly report, gives the two sequences as not identical
      (Relationship ``<>``), since a position on one need not hold on the
      other;
    - any other name, as of a local assembly or a RefSeq contig or scaffold
      (NT_, NW_), is written as it is, and one warning lists these names.

    Returns:
        The name to write for each contig, the problems, and the warnings.
    """
    names: dict[str, str] = {}
    problems: list[str] = []
    local: list[str] = []
    for contig in contigs:
        if REFSEQ_SEQUENCE.match(contig):
            groups = [group for group in aliases or () if contig in group]
            found = sorted({name for group in groups
                            for name in group if GENBANK_SEQUENCE.fullmatch(name)})
            differs = any(contig in getattr(group, "differs", ()) for group in groups)
            if len(found) == 1 and differs:
                source = f"the assembly report {alias_file}" if alias_file else "the alias table"
                problems.append(
                    f"{source} gives the RefSeq sequence '{contig}' and the GenBank sequence "
                    f"'{found[0]}' as not identical (Relationship <>), so a position on "
                    f"'{contig}' need not be the same position on '{found[0]}'. Evaluate and "
                    f"export against the GenBank assembly instead")
                continue
            if len(found) == 1:
                names[contig] = found[0]
                continue
            if aliases is None:
                problems.append(
                    f"the reference sequence '{contig}' is named by its RefSeq accession, "
                    f"while CanSNPer2 2.3.0 names reference sequences by their GenBank "
                    f"accession. Give the names of both with --contig-aliases "
                    f"<assembly_report.txt>, the assembly report NCBI distributes with the "
                    f"assembly")
            elif not found:
                problems.append(
                    f"the alias file {alias_file or 'of --contig-aliases'} gives no GenBank "
                    f"accession for the reference sequence '{contig}'. Add one to it, or use "
                    f"the assembly report NCBI distributes with the assembly "
                    f"(assembly_report.txt), which names both")
            else:
                problems.append(
                    f"the alias file {alias_file or 'of --contig-aliases'} gives more than one "
                    f"GenBank accession for the reference sequence '{contig}': "
                    f"{', '.join(found)}")
            continue
        if not GENBANK_SEQUENCE.fullmatch(contig):
            local.append(contig)
        names[contig] = contig

    written: dict[str, list[str]] = {}
    for contig, name in names.items():
        written.setdefault(name, []).append(contig)
    for name, given in sorted(written.items()):
        if len(given) > 1:
            problems.append(f"the reference sequences {', '.join(given)} would all be written "
                            f"as '{name}'")

    warnings: list[str] = []
    if local:
        other_refseq = ("RefSeq contigs and scaffolds (NT_, NW_) are among them and are not "
                        "replaced. " if any(OTHER_REFSEQ.match(name) for name in local) else "")
        warnings.append(
            f"{len(local)} reference sequences are named neither by a versioned GenBank "
            f"accession nor by a RefSeq accession of the NC_, NZ_ or AC_ kind, which would "
            f"be replaced by its GenBank name; they are written as they are: "
            f"{', '.join(local[:5])}. {other_refseq}CanSNPer2-download will not find a genome "
            f"with sequences of these names, so the reference must be provided by hand")
    return names, problems, warnings


def repeated_positions(scores: Sequence[ScoredSnp]) -> list[str]:
    """Positions that more than one passing record of the scores gives.

    Extract writes one record per alignment column, so a position given twice
    means the scores were joined from more than one run or edited. At two
    nodes it cannot be canonical for both, since one SNP divides the genomes
    in one way only; at one node it would enter the canSNP set twice. Which
    record is right cannot be told from the scores, so this is refused rather
    than resolved by score.
    """
    nodes_at: dict[tuple[str, int], list[str]] = {}
    for score in scores:
        nodes_at.setdefault((score.contig, score.position), []).append(score.node)
    repeated = sorted((key, nodes) for key, nodes in nodes_at.items() if len(nodes) > 1)
    if not repeated:
        return []
    examples = "; ".join(f"{contig}:{position} at {', '.join(nodes)}"
                         for (contig, position), nodes in repeated[:5])
    return [f"the scores give {len(repeated)} positions in more than one passing record, which "
            f"extract does not write: a position marks one branch at most, and only once. "
            f"Were the scores joined from more than one run? {examples}"]


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
    above it that was kept, as that node's children do. Only the genomes of the
    terminal nodes are written, so that a row means one thing, a genome that
    represents that terminal group, as CanSNPer2-data requires; a genome that
    ends up at an inner node is left out with a warning, and it was used to find
    the canSNPs all the same. A genome whose sequence is in the alignment is a
    representative; one placed only through the genome map is a member; an
    aligned genome placed at no node is an outgroup, since extract counted it as
    outside every node. That includes
    the reference genome, the alignment's first record, when it is placed at
    no node; CanSNPer2-data defines no other role for it, so it is listed as
    an outgroup and a warning says so.

    Each genome is written as the assembly accession its name holds. GCA_ and
    GCF_ accessions with the same digits, or two versions of one, are one
    assembly and get one row, whatever order the names come in: the row of a
    representative is kept before that of a member, and a member before an
    outgroup, since the representative's sequence is what the canSNPs were
    found in. One assembly given at more than one node is reported.
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

    warnings: list[str] = []
    no_accession: list[str] = []
    bad_version: list[str] = []
    # Every name given for an assembly, in the order met, as (accession, node, role).
    given: dict[str, list[tuple[str, str, str]]] = {}

    def add(name: str, node: str, role: str) -> None:
        found = ACCESSION.search(name)
        if found is None:
            (bad_version if ANY_ACCESSION.search(name) else no_accession).append(name)
            return
        accession = found.group(0)
        given.setdefault(accession[4:13], []).append((accession, node, role))

    unplaced: list[str] = []
    inner: list[str] = []
    # A node that no kept node gives as its parent is terminal.
    parents_of_kept = {hierarchy.parents.get(node) for node in keeping}
    for genome, node in placed.items():
        if node not in hierarchy:
            unplaced.append(f"{genome} ({node})")
            continue
        at = kept_at(node) or hierarchy.root
        if at in parents_of_kept:
            inner.append(f"{genome} ({at})")
            continue
        add(genome, at, "representative" if genome in matched else "member")
    for record in aligned:
        if record not in aligned_placed:
            add(record, "", "outgroup")

    order = {"representative": 0, "member": 1, "outgroup": 2}
    rows: list[tuple[str, str, str]] = []
    split: list[str] = []
    for entries in given.values():
        rows.append(min(entries, key=lambda entry: order[entry[2]]))
        nodes = sorted({node for _accession, node, _role in entries})
        if len(nodes) > 1:
            names = "/".join(sorted({accession for accession, _node, _role in entries}))
            split.append(f"{names} at {' and '.join(node or 'no node' for node in nodes)}")

    if unplaced:
        warnings.append(
            f"{len(unplaced)} genomes are placed at nodes the hierarchy does not hold and were "
            f"left out of the genomes file: {', '.join(sorted(unplaced)[:5])}")
    if inner:
        warnings.append(
            f"{len(inner)} genomes sit at an inner node and were left out of the genomes file, "
            f"which lists the genomes of the terminal nodes; they were used to find the canSNPs: "
            f"{', '.join(sorted(inner)[:5])}")
    if split:
        warnings.append(
            f"{len(split)} assemblies are given at more than one node, under their GCA_ and "
            f"GCF_ accessions or under two versions; each is listed once, at the node of its "
            f"representative where it has one: {'; '.join(sorted(split)[:5])}")
    if bad_version:
        warnings.append(
            f"{len(bad_version)} genomes hold an assembly accession without a version of one or "
            f"two digits, which CanSNPer2-data requires, and were left out of the genomes file: "
            f"{', '.join(sorted(bad_version)[:5])}")
    if no_accession:
        warnings.append(
            f"{len(no_accession)} genomes hold no assembly accession in their names and were "
            f"left out of the genomes file: {', '.join(sorted(no_accession)[:5])}")
    reference = aligned[0] if aligned else None
    if reference is not None and reference not in aligned_placed:
        found = ACCESSION.search(reference)
        if found is not None and (found.group(0), "", "outgroup") in rows:
            warnings.append(
                f"the reference genome '{reference}' is placed at no node, so the genomes file "
                f"lists it as an outgroup, as extract counted it: outside every node. If it "
                f"belongs to the scheme, place it in the genome map and run extract again")
    rows.sort(key=lambda row: (order[row[2]], row[1], row[0]))
    return rows, warnings


def write_genomes(outdir: str, organism: str, rows: Sequence[tuple[str, str, str]]) -> str:
    path = os.path.join(outdir, organism + GENOMES_SUFFIX)
    with open(path, "w", newline="\n") as handle:
        handle.write("\t".join(GENOME_COLUMNS) + "\n")
        handle.writelines("\t".join(row) + "\n" for row in rows)
    return path


def write_cansnps(outdir: str, organism: str, paths: Sequence[Sequence[str]],
                  chosen: dict[str, Sequence[ScoredSnp]], names: dict[str, str],
                  genome: str, citation: str) -> str:
    """Write every chosen canSNP of each kept node, ranked within the node.

    Each is given by sequence, named as ``names`` gives it, and position on that
    sequence. Rank 1 of a node is its first by :func:`ranking`, and it is the one
    marked typed, which is what CanSNPer2 calls.
    """
    path = os.path.join(outdir, organism + CANSNP_SUFFIX)
    with open(path, "w", newline="\n") as handle:
        handle.write("\t".join(CANSNP_COLUMNS) + "\n")
        for tree_path in paths:
            node = tree_path[-1]
            handle.writelines("\t".join([
                node, str(rank), citation, genome, names[score.contig], str(score.position),
                score.derived, score.ancestral, str(score.total), "yes" if rank == 1 else "no",
            ]) + "\n" for rank, score in enumerate(chosen[node], start=1))
    return path


def write_scheme(outdir: str, organism: str, paths: Sequence[Sequence[str]],
                 best: dict[str, ScoredSnp], names: dict[str, str], genome: str, strain: str,
                 citation: str, genbank_id: str, refseq_id: str,
                 assembly_name: str) -> list[str]:
    """Write the tree and references files and return their paths.

    The canSNPs are written by :func:`write_cansnps`. ``names`` gives the name each
    contig of the reference is written under.
    """
    os.makedirs(outdir, exist_ok=True)
    written = []

    tree_path = os.path.join(outdir, organism + TREE_SUFFIX)
    with open(tree_path, "w", newline="\n") as handle:
        handle.writelines("\t".join(path) + "\n" for path in paths)
    written.append(tree_path)

    references_path = os.path.join(outdir, organism + REFERENCES_SUFFIX)
    with open(references_path, "w", newline="\n") as handle:
        handle.write("\t".join(REFERENCE_COLUMNS) + "\n")
        handle.write("\t".join([genome, strain, genbank_id, refseq_id, assembly_name]) + "\n")
    written.append(references_path)

    return written


def main(scores=None, hierarchy=None, reference=None, outdir=".", organism=None,
         genome=None, strain=None, citation=None, genbank_id=None, refseq_id=None,
         assembly_name=None, prune=False, genome_map=None, alignment=None,
         snps_per_node=1, spacing=1000, contig_aliases=None) -> int:
    """Write a scored scheme as CanSNPer2-data source files. Returns an exit status.

    ``snps_per_node`` is how many canSNPs of each node to keep, at least
    ``spacing`` bases apart, or None for all of them. With more than one the
    full set is written to <organism>_canSNPs.txt. ``contig_aliases`` is an
    NCBI assembly report or UCSC chromAlias table that gives the GenBank
    accession of each sequence the reference names by its RefSeq accession.
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
        aliases = read_contig_aliases(contig_aliases) if contig_aliases else None
    except (OSError, ValueError, HierarchyError, ReferenceError) as error:
        logging.error(error)
        return 1

    chosen = select_per_node(scored, snps_per_node, spacing)
    best = {node: picked[0] for node, picked in chosen.items() if picked}
    problems, warnings = check_scheme(tree, best, sequence)
    problems += repeated_positions(scored)
    # Only the sequences that carry a canSNP of the scheme are named in it.
    used = {score.contig for picked in chosen.values() for score in picked}
    names, naming_problems, naming_warnings = sequence_names(
        [name for name in sequence.names if name in used], aliases, contig_aliases)
    problems += naming_problems
    warnings += naming_warnings
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
            message = (
                f"{len(without)} nodes have no canSNP and were left out of the scheme, so their "
                f"children now hang from the nearest node that has one: {', '.join(without[:5])}")
            if clusters is not None:
                moved = sum(len(clusters.get(node, [])) for node in without)
                depth = max(len(tree.path(node)) for node in without)
                message += (f". {moved} genomes of the genome map sat at these nodes and now sit "
                            f"at an ancestor; the deepest pruned depth is {depth} (the root is "
                            f"depth 1)")
            warnings.append(message)
        else:
            problems.append(
                f"{len(without)} nodes have no canSNP: {', '.join(without[:5])}. Every node of a "
                f"CanSNPer2 scheme carries exactly one, so either lower the evaluate cutoff or pass "
                f"--prune to leave these nodes out")

    if not report(problems, warnings, "scheme"):
        return 2

    kept = [node for node in tree.nodes if node in best]
    paths = paths_for(tree, kept)
    written = write_scheme(outdir, organism, paths, best, names, genome, strain,
                           citation, genbank_id, refseq_id, assembly_name)
    written.append(write_cansnps(outdir, organism, paths, chosen, names, genome, citation))
    # The optional files of an earlier export into the same directory would
    # otherwise stay beside this scheme, and describe another one. The SNP file of
    # CanSNPer2 2.3.0 and earlier is no longer written at all.
    for suffix, wanted in ((SNP_SUFFIX, False), (GENOMES_SUFFIX, clusters is not None)):
        stale = os.path.join(outdir, organism + suffix)
        if not wanted and os.path.exists(stale):
            os.remove(stale)
            logging.warning(f"Removed {stale}, which an earlier export wrote and this one does "
                            f"not")
    if snps_per_node != 1:
        counts = sorted(len(chosen[path[-1]]) for path in paths)
        short = [path[-1] for path in paths
                 if snps_per_node is not None and len(chosen[path[-1]]) < snps_per_node]
        logging.info(
            f"The canSNP set holds {sum(counts)} canSNPs, {counts[0]} to {counts[-1]} per node")
        shared = [path[-1] for path in paths
                  if len({(score.contig, score.locus) for score in chosen[path[-1]]
                          if score.locus})
                  < sum(1 for score in chosen[path[-1]] if score.locus)]
        if shared:
            logging.warning(
                f"{len(shared)} nodes hold canSNPs that share a gene, which one deletion or "
                f"recombinant segment would take together: {', '.join(shared[:5])}")
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
