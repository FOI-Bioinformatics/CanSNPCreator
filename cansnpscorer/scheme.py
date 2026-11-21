"""Reading a published CanSNPer2-data scheme back in.

The export subcommand writes the three source files of a scheme; this module
reads them again, so that an existing scheme can be examined without being
rebuilt. The reader is deliberately tolerant: a scheme is read here in order to
be audited, so what is wrong with it must survive reading rather than stop it.
Rows that cannot be read at all are collected as problems and handed to the
caller together with the scheme.

The three files are described in the schemes/README.md of the CanSNPer2-data
repository. In short: the tree file holds one line per path from the root, the
SNP file one line per node with a header, and the references file one line per
reference genome with a header.

The SNP file and the optional canSNP set are read in the format of CanSNPer2
2.3.0, in which each canSNP names a reference sequence, by its GenBank
accession, and gives a 1-based position on that sequence. Columns are found by
the names in the header. A file written before 2.3.0, with one position along
the reference sequences joined in file order, is refused with the command that
converts it, since that position cannot be placed on a sequence without the
reference genome.
"""

import glob
import logging
import os
from collections.abc import Sequence
from dataclasses import dataclass, field

from cansnpscorer.hierarchy import Hierarchy, HierarchyError

TREE_SUFFIX = "_tree.txt"
SNP_SUFFIX = "_snp.txt"
REFERENCES_SUFFIX = "_references.txt"
GENOMES_SUFFIX = "_genomes.txt"
SNP_SET_SUFFIX = "_canSNPs.txt"
SNP_SET_COLUMNS = ("node", "rank", "genome", "sequence", "sequence_position", "derived_base",
                   "ancestral_base", "score")

SNP_COLUMNS = ("snp_id", "strain", "reference", "genome", "sequence", "sequence_position",
               "derived_base", "ancestral_base")
#: What a file with the header of CanSNPer2 before 2.3.0 is told to do.
CONVERT_HINT = ("it has the header of CanSNPer2 before 2.3.0, with one position along the "
                "reference sequences joined in file order; convert it with "
                "CanSNPer2-database --convert (CanSNPer2 2.3.0), which gives each canSNP a "
                "sequence and a position on it")
REFERENCE_COLUMNS = ("genome", "strain", "genbank_id", "refseq_id", "assembly_name")
GENOME_COLUMNS = ("genome", "node", "role")
GENOME_ROLES = ("representative", "member", "outgroup")


class SchemeError(ValueError):
    """The scheme's files are not there, or cannot be read at all."""


@dataclass(frozen=True)
class SchemeSnp:
    """One row of the SNP file."""

    node: str
    strain: str
    reference: str
    genome: str
    sequence: str
    """The reference sequence the canSNP is on, as the SNP file names it."""
    position: int
    """1-based position on ``sequence``."""
    derived: str
    ancestral: str
    line: int = 0
    rank: int = 1
    """Its rank among the canSNPs of its node; the SNP file holds rank 1."""

    @property
    def where(self) -> str:
        return f"{self.genome} {self.sequence}:{self.position}"

    def tag(self, text: str) -> str:
        """A finding's text, saying which canSNP of the node it is about when
        it is not the one the SNP file gives."""
        return text if self.rank == 1 else f"the canSNP of rank {self.rank}: {text}"


@dataclass(frozen=True)
class SchemeReference:
    """One row of the references file."""

    genome: str
    strain: str
    genbank_id: str
    refseq_id: str
    assembly_name: str
    line: int = 0


@dataclass
class Scheme:
    """A scheme as its three source files hold it."""

    organism: str
    paths: list[list[str]] = field(default_factory=list)
    snps: list[SchemeSnp] = field(default_factory=list)
    references: list[SchemeReference] = field(default_factory=list)
    hierarchy: Hierarchy | None = None
    genomes: list[tuple[str, str, str]] = field(default_factory=list)
    """The genomes the scheme was made from, as (accession, node, role), when
    the scheme lists them in its optional genomes file."""
    snp_set: list[SchemeSnp] = field(default_factory=list)
    """Every canSNP of each node, ranked, when the scheme lists more than one
    in its optional canSNP set; rank 1 is the SNP of the SNP file."""
    # Nodes placed under more than one parent, as (node, parents). The tree
    # reader of CanSNPer2 does not refuse these; it creates a second node with
    # the same name and attaches the SNP to one of them.
    conflicts: list[tuple[str, list[str]]] = field(default_factory=list)
    problems: list[str] = field(default_factory=list)

    @property
    def nodes(self) -> list[str]:
        """Every node of the tree, parents before their children where known."""
        if self.hierarchy is not None:
            return self.hierarchy.nodes
        seen: list[str] = []
        for path in self.paths:
            for node in path:
                if node not in seen:
                    seen.append(node)
        return seen

    @property
    def snps_by_node(self) -> dict[str, list[SchemeSnp]]:
        by_node: dict[str, list[SchemeSnp]] = {}
        for snp in self.snps:
            by_node.setdefault(snp.node, []).append(snp)
        return by_node

    @property
    def checked_snps(self) -> list[SchemeSnp]:
        """The SNPs to check against a genome or an alignment: those of the SNP
        file, and every further rank of the canSNP set."""
        return self.snps + [snp for snp in self.snp_set if snp.rank > 1]

    @property
    def references_by_genome(self) -> dict[str, SchemeReference]:
        return {reference.genome: reference for reference in self.references}


def schemes_in(directory: str) -> list[str]:
    """The organisms whose tree file is in this directory."""
    found = glob.glob(os.path.join(glob.escape(directory), "*" + TREE_SUFFIX))
    return sorted(os.path.basename(path)[:-len(TREE_SUFFIX)] for path in found)


def paths_for(directory: str, organism: str) -> tuple[str, str, str]:
    """The three file paths of one scheme."""
    return (
        os.path.join(directory, organism + TREE_SUFFIX),
        os.path.join(directory, organism + SNP_SUFFIX),
        os.path.join(directory, organism + REFERENCES_SUFFIX),
    )


def _rows(path: str, columns: Sequence[str], problems: list[str]) -> list[tuple[int, list[str]]]:
    """Read a tab separated file with a header, by the names of its columns.

    Each row is returned with its fields in the order of ``columns``, whatever
    order the file gives them in. A header without one of the columns is a
    problem and no rows are read; a column not among them is a problem, and the
    rows are read without it. A header with a position and no
    sequence_position is that of CanSNPer2 before 2.3.0, which is refused.

    Raises:
        SchemeError: If the file has the header of CanSNPer2 before 2.3.0.
    """
    rows: list[tuple[int, list[str]]] = []
    with open(path, newline="") as handle:
        # Only LF and CR LF end a line; splitlines would also split on Unicode
        # line separators inside a field.
        lines = [line.rstrip("\r") for line in handle.read().split("\n")]
        while lines and not lines[-1]:
            lines.pop()

    if not lines:
        problems.append(f"{os.path.basename(path)} is empty")
        return rows

    name = os.path.basename(path)
    header = [column.strip() for column in lines[0].split("\t")]
    if "sequence_position" in columns and "position" in header \
            and "sequence_position" not in header:
        raise SchemeError(f"{name}: {CONVERT_HINT}")
    missing = [column for column in columns if column not in header]
    if missing or len(set(header)) != len(header):
        problems.append(
            f"{name} has the header {', '.join(header)} where the columns "
            f"{', '.join(columns)} are expected, each once")
        return rows
    extra = [column for column in header if column not in columns]
    if extra:
        problems.append(f"{name} has the columns {', '.join(extra)}, which are not read")
    order = [header.index(column) for column in columns]

    for number, line in enumerate(lines[1:], 2):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != len(header):
            problems.append(
                f"{name} line {number} holds {len(fields)} columns where "
                f"{len(header)} are expected")
            continue
        rows.append((number, [fields[index] for index in order]))
    return rows


def read_tree_file(path: str) -> tuple[list[list[str]], list[str]]:
    """Read the tree file as a list of paths from the root."""
    problems: list[str] = []
    paths: list[list[str]] = []

    with open(path, newline="") as handle:
        for number, line in enumerate(
                (line.rstrip("\r") for line in handle.read().split("\n")), 1):
            if not line.strip():
                continue
            nodes = [field.strip() for field in line.split("\t")]
            if any(not node for node in nodes):
                problems.append(
                    f"{os.path.basename(path)} line {number} holds an empty node name")
                continue
            paths.append(nodes)

    if not paths:
        problems.append(f"{os.path.basename(path)} names no nodes")
    return paths, problems


def hierarchy_from_paths(paths: Sequence[Sequence[str]]) -> tuple[Hierarchy | None,
                                                                 list[tuple[str, list[str]]],
                                                                 list[str]]:
    """Turn the tree file's paths into a hierarchy.

    A node named under two parents is reported rather than raised on, and the
    first parent is kept, so that the rest of the scheme can still be examined.
    """
    parents: dict[str, str | None] = {}
    seen: dict[str, list[str]] = {}
    problems: list[str] = []

    for path in paths:
        for depth, node in enumerate(path):
            parent = path[depth - 1] if depth else None
            label = parent if parent is not None else "the root"
            if label not in seen.setdefault(node, []):
                seen[node].append(label)
            parents.setdefault(node, parent)

    conflicts = [(node, found) for node, found in seen.items() if len(found) > 1]

    try:
        hierarchy = Hierarchy(list(parents.items()))
    except HierarchyError as error:
        problems.append(f"the tree file is not a single rooted tree: {error}")
        hierarchy = None

    return hierarchy, sorted(conflicts), problems


def read_scheme(directory: str, organism: str | None = None) -> Scheme:
    """Read one scheme from a directory of CanSNPer2-data source files.

    Args:
        directory: The directory holding the three files.
        organism: Which scheme to read. With one scheme in the directory it may
            be left out.

    Returns:
        The scheme, with whatever could not be read listed in ``problems``.

    Raises:
        SchemeError: If the scheme's files are not there, or its SNP file or
            canSNP set has the header of CanSNPer2 before 2.3.0.
        OSError: If a file cannot be read.
    """
    available = schemes_in(directory)
    if organism is None:
        if len(available) == 1:
            organism = available[0]
        elif not available:
            raise SchemeError(f"{directory} holds no *{TREE_SUFFIX} file")
        else:
            raise SchemeError(
                f"{directory} holds {len(available)} schemes, so one must be named with "
                f"--organism: {', '.join(available)}")

    tree_path, snp_path, references_path = paths_for(directory, organism)
    missing = [path for path in (tree_path, snp_path, references_path)
               if not os.path.isfile(path)]
    if missing:
        raise SchemeError(
            f"these files of the '{organism}' scheme are missing: "
            f"{', '.join(os.path.basename(path) for path in missing)}")

    logging.info(f"Reading the '{organism}' scheme from {directory}")

    paths, problems = read_tree_file(tree_path)
    hierarchy, conflicts, tree_problems = hierarchy_from_paths(paths)
    problems += tree_problems

    snps: list[SchemeSnp] = []
    for number, fields in _rows(snp_path, SNP_COLUMNS, problems):
        node, strain, reference, genome, sequence, position, derived, ancestral = fields
        try:
            value = int(position)
        except ValueError:
            problems.append(
                f"{os.path.basename(snp_path)} line {number} gives the sequence_position of "
                f"'{node}' as '{position}', which is not a whole number")
            continue
        if not sequence.strip():
            problems.append(
                f"{os.path.basename(snp_path)} line {number} names no sequence for '{node}'")
            continue
        snps.append(SchemeSnp(node=node.strip(), strain=strain, reference=reference,
                              genome=genome.strip(), sequence=sequence.strip(), position=value,
                              derived=derived.strip(), ancestral=ancestral.strip(),
                              line=number))

    references: list[SchemeReference] = []
    for number, fields in _rows(references_path, REFERENCE_COLUMNS, problems):
        references.append(SchemeReference(*[field.strip() for field in fields], line=number))

    genomes: list[tuple[str, str, str]] = []
    genomes_path = os.path.join(directory, organism + GENOMES_SUFFIX)
    if os.path.isfile(genomes_path):
        for number, fields in _rows(genomes_path, GENOME_COLUMNS, problems):
            genome, node, role = (field.strip() for field in fields)
            if role not in GENOME_ROLES:
                problems.append(f"{os.path.basename(genomes_path)} line {number} gives the "
                                f"role '{role}'")
                continue
            genomes.append((genome, node, role))

    snp_set: list[SchemeSnp] = []
    snp_set_path = os.path.join(directory, organism + SNP_SET_SUFFIX)
    if os.path.isfile(snp_set_path):
        snp_set = read_snp_set(snp_set_path, snps, problems)

    logging.info(
        f"Read {len(paths)} tree lines, {len(snps)} SNPs and {len(references)} reference genomes"
        + (f", {len(genomes)} genomes the scheme was made from" if genomes else "")
        + (f", and a set of {len(snp_set)} ranked canSNPs" if snp_set else ""))

    return Scheme(organism=organism, paths=paths, snps=snps, references=references,
                  hierarchy=hierarchy, genomes=genomes, snp_set=snp_set, conflicts=conflicts,
                  problems=problems)


def read_snp_set(path: str, snps: Sequence[SchemeSnp], problems: list[str]) -> list[SchemeSnp]:
    """Read the optional ranked canSNPs of each node.

    Rank 1 of a node must be the SNP the SNP file gives it, since that is the
    one CanSNPer2 calls, and the ranks of a node must run 1, 2, 3 ... without
    a gap or a repeat; a deviation is a problem of the scheme.
    """
    name = os.path.basename(path)
    snp_set: list[SchemeSnp] = []
    for number, fields in _rows(path, SNP_SET_COLUMNS, problems):
        node, rank, genome, sequence, position, derived, ancestral, _score = (
            f.strip() for f in fields)
        try:
            rank_value, position_value = int(rank), int(position)
        except ValueError:
            problems.append(f"{name} line {number} gives the rank '{rank}' and sequence_position "
                            f"'{position}', which are not both whole numbers")
            continue
        if not sequence:
            problems.append(f"{name} line {number} names no sequence for '{node}'")
            continue
        snp_set.append(SchemeSnp(node=node, strain="", reference="", genome=genome,
                                 sequence=sequence, position=position_value, derived=derived,
                                 ancestral=ancestral, line=number, rank=rank_value))

    # Export numbers each node's canSNPs 1, 2, 3 ... in file order. A gap or a
    # repeat means the file was edited or joined, and a repeated rank 1 leaves
    # open which canSNP the SNP file's is meant to be.
    ranks: dict[str, list[int]] = {}
    for snp in snp_set:
        ranks.setdefault(snp.node, []).append(snp.rank)
    for node, found in sorted(ranks.items()):
        if found != list(range(1, len(found) + 1)):
            problems.append(
                f"{name} gives node '{node}' the ranks {', '.join(map(str, found))}, where "
                f"1 to {len(found)} in order are expected")

    called = {snp.node: (snp.genome, snp.sequence, snp.position, snp.derived, snp.ancestral)
              for snp in snps}
    first = {snp.node: (snp.genome, snp.sequence, snp.position, snp.derived, snp.ancestral)
             for snp in snp_set if snp.rank == 1}
    for node in sorted(set(called) & set(first)):
        if called[node] != first[node]:
            problems.append(f"{name} gives node '{node}' a rank 1 that is not the SNP the SNP "
                            f"file gives it")
    for node in sorted(set(called) - set(first)):
        problems.append(f"{name} has no canSNP of rank 1 for node '{node}'")
    for node in sorted(set(first) - set(called)):
        problems.append(f"{name} ranks canSNPs for node '{node}', which the SNP file does "
                        f"not hold")
    return snp_set
