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
CANSNP_SUFFIX = "_canSNPs.txt"
#: The canSNP file of a scheme: every canSNP of every node, ranked within the node, with
#: ``typed`` marking the ones CanSNPer2 calls. ``score`` is empty for a canSNP taken from
#: a publication.
CANSNP_COLUMNS = ("node", "rank", "reference", "genome", "sequence", "sequence_position",
                  "derived_base", "ancestral_base", "score", "typed")
#: The SNP file of CanSNPer2 2.3.0, one canSNP per node, which a scheme exported before the
#: canSNP file still has. It is read as rank 1 of its node, typed and without a score.
SNP_COLUMNS = ("snp_id", "strain", "reference", "genome", "sequence", "sequence_position",
               "derived_base", "ancestral_base")
TYPED_VALUES = {"yes": True, "no": False}
#: What a file with the header of CanSNPer2 before 2.3.0 is told to do.
CONVERT_HINT = ("it has the header of CanSNPer2 before 2.3.0, with one position along the "
                "reference sequences joined in file order; convert it with "
                "CanSNPer2-database --convert (CanSNPer2 2.3.0), which gives each canSNP a "
                "sequence and a position on it")
REFERENCE_COLUMNS = ("genome", "strain", "genbank_id", "refseq_id", "assembly_name")
GENOME_COLUMNS = ("genome", "node", "role")
#: The genomes file may record where each genome comes from, in two more columns.
GENOME_COLUMNS_WITH_SOURCE = GENOME_COLUMNS + ("basis", "source")
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
    """Its rank among the canSNPs of its node; rank 1 is the node's best canSNP."""
    score: str = ""
    """What it scored when the scheme was made; empty for a published canSNP."""
    typed: bool = True
    """Whether the scheme types it, which is what CanSNPer2 calls."""

    @property
    def where(self) -> str:
        return f"{self.genome} {self.sequence}:{self.position}"

    def tag(self, text: str) -> str:
        """A finding's text, saying which canSNP of the node it is about when
        it is not the one the scheme types."""
        return text if self.typed else f"the canSNP of rank {self.rank}: {text}"


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
    """Every canSNP of each node, ranked within the node, as the canSNP file
    holds them. ``snps`` is the subset the scheme types."""
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
        """The canSNPs to check against a genome or an alignment: every canSNP of
        every node, the alternatives included."""
        return list(self.snp_set) if self.snp_set else list(self.snps)

    @property
    def references_by_genome(self) -> dict[str, SchemeReference]:
        return {reference.genome: reference for reference in self.references}


def schemes_in(directory: str) -> list[str]:
    """The organisms whose tree file is in this directory."""
    found = glob.glob(os.path.join(glob.escape(directory), "*" + TREE_SUFFIX))
    return sorted(os.path.basename(path)[:-len(TREE_SUFFIX)] for path in found)


def paths_for(directory: str, organism: str) -> tuple[str, str, str]:
    """The three file paths of one scheme: its tree, its canSNPs and its references.

    A scheme exported before the canSNP file keeps its canSNPs in ``_snp.txt``, one
    per node; that file is used when there is no ``_canSNPs.txt``.
    """
    cansnps = os.path.join(directory, organism + CANSNP_SUFFIX)
    if not os.path.isfile(cansnps):
        old = os.path.join(directory, organism + SNP_SUFFIX)
        if os.path.isfile(old):
            cansnps = old
    return (
        os.path.join(directory, organism + TREE_SUFFIX),
        cansnps,
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

    snp_set = read_cansnps(snp_path, problems)
    snps = [snp for snp in snp_set if snp.typed]

    references: list[SchemeReference] = []
    for number, fields in _rows(references_path, REFERENCE_COLUMNS, problems):
        references.append(SchemeReference(*[field.strip() for field in fields], line=number))

    genomes: list[tuple[str, str, str]] = []
    genomes_path = os.path.join(directory, organism + GENOMES_SUFFIX)
    if os.path.isfile(genomes_path):
        # basis and source are optional, and both or neither; they say where a genome
        # comes from and the audit does not read them.
        with open(genomes_path, newline="") as handle:
            first = handle.readline().rstrip("\r\n").split("\t")
        columns = (GENOME_COLUMNS_WITH_SOURCE if "basis" in first or "source" in first
                   else GENOME_COLUMNS)
        for number, fields in _rows(genomes_path, columns, problems):
            genome, node, role = (field.strip() for field in fields[:3])
            if role not in GENOME_ROLES:
                problems.append(f"{os.path.basename(genomes_path)} line {number} gives the "
                                f"role '{role}'")
                continue
            genomes.append((genome, node, role))

    logging.info(
        f"Read {len(paths)} tree lines, {len(snps)} SNPs and {len(references)} reference genomes"
        + (f", {len(genomes)} genomes the scheme was made from" if genomes else "")
        + (f", and a set of {len(snp_set)} ranked canSNPs" if snp_set else ""))

    return Scheme(organism=organism, paths=paths, snps=snps, references=references,
                  hierarchy=hierarchy, genomes=genomes, snp_set=snp_set, conflicts=conflicts,
                  problems=problems)


def read_cansnps(path: str, problems: list[str]) -> list[SchemeSnp]:
    """Read the canSNPs of a scheme: every canSNP of every node, ranked within the node.

    A file with the header of CanSNPer2 2.3.0, one canSNP per node, is read as rank 1 of
    its node, typed and without a score. The ranks of a node must run 1, 2, 3 ... without
    a gap or a repeat, and exactly one canSNP of a node is typed, since CanSNPer2 calls
    one canSNP per node; a deviation is a problem of the scheme.
    """
    name = os.path.basename(path)
    old_format = SNP_SUFFIX if path.endswith(SNP_SUFFIX) else None
    columns = SNP_COLUMNS if old_format else CANSNP_COLUMNS
    cansnps: list[SchemeSnp] = []
    for number, fields in _rows(path, columns, problems):
        values = [field.strip() for field in fields]
        if old_format:
            node, _strain, reference, genome, sequence, position, derived, ancestral = values
            rank, score, typed = "1", "", "yes"
        else:
            (node, rank, reference, genome, sequence, position, derived, ancestral,
             score, typed) = values
        try:
            rank_value, position_value = int(rank), int(position)
        except ValueError:
            problems.append(f"{name} line {number} gives the rank '{rank}' and sequence_position "
                            f"'{position}', which are not both whole numbers")
            continue
        if not sequence:
            problems.append(f"{name} line {number} names no sequence for '{node}'")
            continue
        if typed.lower() not in TYPED_VALUES:
            problems.append(f"{name} line {number} gives typed as '{typed}' for '{node}', "
                            f"where yes or no is expected")
            continue
        cansnps.append(SchemeSnp(node=node, strain="", reference=reference, genome=genome,
                                 sequence=sequence, position=position_value, derived=derived,
                                 ancestral=ancestral, line=number, rank=rank_value,
                                 score=score, typed=TYPED_VALUES[typed.lower()]))

    # Export numbers each node's canSNPs 1, 2, 3 ... in file order. A gap or a
    # repeat means the file was edited or joined.
    ranks: dict[str, list[int]] = {}
    typed_count: dict[str, int] = {}
    for snp in cansnps:
        ranks.setdefault(snp.node, []).append(snp.rank)
        typed_count[snp.node] = typed_count.get(snp.node, 0) + (1 if snp.typed else 0)
    for node, found in sorted(ranks.items()):
        if found != list(range(1, len(found) + 1)):
            problems.append(
                f"{name} gives node '{node}' the ranks {', '.join(map(str, found))}, where "
                f"1 to {len(found)} in order are expected")
    for node, count in sorted(typed_count.items()):
        if count == 0:
            problems.append(f"{name} marks no canSNP of node '{node}' typed, so CanSNPer2 "
                            f"would not call it")
        elif count > 1:
            problems.append(f"{name} marks {count} canSNPs of node '{node}' typed, where "
                            f"CanSNPer2 calls one canSNP per node")
    return cansnps
