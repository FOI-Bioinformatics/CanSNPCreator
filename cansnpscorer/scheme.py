"""Reading a published CanSNPer2-data scheme back in.

The export subcommand writes the three source files of a scheme; this module
reads them again, so that an existing scheme can be examined without being
rebuilt. The reader is deliberately tolerant: a scheme is read here in order to
be audited, so what is wrong with it must survive reading rather than stop it.
Rows that cannot be read at all are collected as problems and handed to the
caller together with the scheme.

The three files are described in the source/README.md of the CanSNPer2-data
repository. In short: the tree file holds one line per path from the root, the
SNP file one line per node with a header, and the references file one line per
reference genome with a header.
"""

import glob
import logging
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from cansnpscorer.hierarchy import Hierarchy, HierarchyError

TREE_SUFFIX = "_tree.txt"
SNP_SUFFIX = "_snp.txt"
REFERENCES_SUFFIX = "_references.txt"
GENOMES_SUFFIX = "_genomes.txt"

SNP_COLUMNS = ("snp_id", "strain", "reference", "genome", "position",
               "derived_base", "ancestral_base")
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
    position: int
    derived: str
    ancestral: str
    line: int = 0

    @property
    def where(self) -> str:
        return f"{self.genome}:{self.position}"


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
    paths: List[List[str]] = field(default_factory=list)
    snps: List[SchemeSnp] = field(default_factory=list)
    references: List[SchemeReference] = field(default_factory=list)
    hierarchy: Optional[Hierarchy] = None
    genomes: List[Tuple[str, str, str]] = field(default_factory=list)
    """The genomes the scheme was made from, as (accession, node, role), when
    the scheme lists them in its optional genomes file."""
    # Nodes placed under more than one parent, as (node, parents). The tree
    # reader of CanSNPer2 does not refuse these; it creates a second node with
    # the same name and attaches the SNP to one of them.
    conflicts: List[Tuple[str, List[str]]] = field(default_factory=list)
    problems: List[str] = field(default_factory=list)

    @property
    def nodes(self) -> List[str]:
        """Every node of the tree, parents before their children where known."""
        if self.hierarchy is not None:
            return self.hierarchy.nodes
        seen: List[str] = []
        for path in self.paths:
            for node in path:
                if node not in seen:
                    seen.append(node)
        return seen

    @property
    def snps_by_node(self) -> Dict[str, List[SchemeSnp]]:
        by_node: Dict[str, List[SchemeSnp]] = {}
        for snp in self.snps:
            by_node.setdefault(snp.node, []).append(snp)
        return by_node

    @property
    def references_by_genome(self) -> Dict[str, SchemeReference]:
        return {reference.genome: reference for reference in self.references}


def schemes_in(directory: str) -> List[str]:
    """The organisms whose tree file is in this directory."""
    found = glob.glob(os.path.join(directory, "*" + TREE_SUFFIX))
    return sorted(os.path.basename(path)[:-len(TREE_SUFFIX)] for path in found)


def paths_for(directory: str, organism: str) -> Tuple[str, str, str]:
    """The three file paths of one scheme."""
    return (
        os.path.join(directory, organism + TREE_SUFFIX),
        os.path.join(directory, organism + SNP_SUFFIX),
        os.path.join(directory, organism + REFERENCES_SUFFIX),
    )


def _rows(path: str, columns: Sequence[str], problems: List[str]) -> List[Tuple[int, List[str]]]:
    """Read a tab separated file with a header, checking the column names."""
    rows: List[Tuple[int, List[str]]] = []
    with open(path) as handle:
        lines = handle.read().splitlines()

    if not lines:
        problems.append(f"{os.path.basename(path)} is empty")
        return rows

    header = lines[0].split("\t")
    if [name.strip() for name in header] != list(columns):
        problems.append(
            f"{os.path.basename(path)} has the header {', '.join(header)} where "
            f"{', '.join(columns)} is expected")
        return rows

    for number, line in enumerate(lines[1:], 2):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != len(columns):
            problems.append(
                f"{os.path.basename(path)} line {number} holds {len(fields)} columns where "
                f"{len(columns)} are expected")
            continue
        rows.append((number, fields))
    return rows


def read_tree_file(path: str) -> Tuple[List[List[str]], List[str]]:
    """Read the tree file as a list of paths from the root."""
    problems: List[str] = []
    paths: List[List[str]] = []

    with open(path) as handle:
        for number, line in enumerate(handle.read().splitlines(), 1):
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


def hierarchy_from_paths(paths: Sequence[Sequence[str]]) -> Tuple[Optional[Hierarchy],
                                                                 List[Tuple[str, List[str]]],
                                                                 List[str]]:
    """Turn the tree file's paths into a hierarchy.

    A node named under two parents is reported rather than raised on, and the
    first parent is kept, so that the rest of the scheme can still be examined.
    """
    parents: Dict[str, Optional[str]] = {}
    seen: Dict[str, List[str]] = {}
    problems: List[str] = []

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


def read_scheme(directory: str, organism: Optional[str] = None) -> Scheme:
    """Read one scheme from a directory of CanSNPer2-data source files.

    Args:
        directory: The directory holding the three files.
        organism: Which scheme to read. With one scheme in the directory it may
            be left out.

    Returns:
        The scheme, with whatever could not be read listed in ``problems``.

    Raises:
        SchemeError: If the scheme's files are not there.
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

    snps: List[SchemeSnp] = []
    for number, fields in _rows(snp_path, SNP_COLUMNS, problems):
        node, strain, reference, genome, position, derived, ancestral = fields
        try:
            value = int(position)
        except ValueError:
            problems.append(
                f"{os.path.basename(snp_path)} line {number} gives the position of '{node}' as "
                f"'{position}', which is not a whole number")
            continue
        snps.append(SchemeSnp(node=node.strip(), strain=strain, reference=reference,
                              genome=genome.strip(), position=value,
                              derived=derived.strip(), ancestral=ancestral.strip(),
                              line=number))

    references: List[SchemeReference] = []
    for number, fields in _rows(references_path, REFERENCE_COLUMNS, problems):
        references.append(SchemeReference(*[field.strip() for field in fields], line=number))

    genomes: List[Tuple[str, str, str]] = []
    genomes_path = os.path.join(directory, organism + GENOMES_SUFFIX)
    if os.path.isfile(genomes_path):
        for number, fields in _rows(genomes_path, GENOME_COLUMNS, problems):
            genome, node, role = (field.strip() for field in fields)
            if role not in GENOME_ROLES:
                problems.append(f"{os.path.basename(genomes_path)} line {number} gives the "
                                f"role '{role}'")
                continue
            genomes.append((genome, node, role))

    logging.info(
        f"Read {len(paths)} tree lines, {len(snps)} SNPs and {len(references)} reference genomes"
        + (f", and {len(genomes)} genomes the scheme was made from" if genomes else ""))

    return Scheme(organism=organism, paths=paths, snps=snps, references=references,
                  hierarchy=hierarchy, genomes=genomes, conflicts=conflicts, problems=problems)
