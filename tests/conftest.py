"""Shared test data for CanSNPCreator.

The ``dataset`` fixture writes a small but complete set of input files and
records the canSNPs that were planted in them, so a test can assert what the
tools are expected to recover.

The hierarchy is::

    root
    |-- A
    |   |-- A.1   (s1, s1b)
    |   `-- A.2   (s2)
    `-- B
        |-- B.1   (s3)
        `-- B.2   (s4)

The alignment holds the reference first, as the tools require, and covers a
two-contig genome: columns 0-699 are contig ``c1`` and columns 700-999 are
contig ``c2``.
"""

import os
import random
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

GENOMES = ["ref", "s1", "s1b", "s2", "s3", "s4"]
CLUSTERS = {"A.1": ["s1", "s1b"], "A.2": ["s2"], "B.1": ["s3"], "B.2": ["s4"]}
HIERARCHY = [{"root": [{"A": ["A.1", "A.2"]}, {"B": ["B.1", "B.2"]}]}]

ALIGNMENT_LENGTH = 1000
CONTIGS = [("c1", 0, 700), ("c2", 700, 1000)]

# node -> (alignment column, ancestral base, derived base)
PLANTED: Dict[str, Tuple[int, str, str]] = {
    "root": (100, "A", "T"),
    "A": (200, "G", "C"),
    "B": (300, "T", "A"),
    "A.1": (400, "C", "G"),
    "A.2": (500, "G", "A"),
    "B.1": (600, "T", "C"),
    "B.2": (700, "A", "T"),
}

# A SNP shared by two isolates on different branches: a SNP, but not canonical.
INCONGRUENT = (800, "A", "G", ["s1", "s3"])
# Sites the extractor is expected to leave out.
TRIALLELIC = 850
GAPPED = 900

# Genes as they appear in a GFF: 1-based, local to the contig.
GENES = [
    ("c1", 151, 450, "+"),
    ("c1", 551, 650, "-"),
    ("c2", 11, 200, "+"),
]
# Tandem repeats as TRF -ngs reports them: 1-based, local to the contig.
REPEATS = [("c1", 60, 95), ("c2", 120, 160)]


def node_isolates(node: str) -> List[str]:
    """The isolates below a node, which is what a canSNP for it must separate."""
    descendants = {
        "root": ["A.1", "A.2", "B.1", "B.2"],
        "A": ["A.1", "A.2"],
        "B": ["B.1", "B.2"],
    }
    clusters = descendants.get(node, [node])
    return [isolate for cluster in clusters for isolate in CLUSTERS[cluster]]


def build_alignment() -> Dict[str, str]:
    """Build the alignment, planting one canSNP per node."""
    random.seed(20261112)
    backbone = [random.choice("ACGT") for _ in range(ALIGNMENT_LENGTH)]
    sequences = {name: list(backbone) for name in GENOMES}

    for node, (column, ancestral, derived) in PLANTED.items():
        for name in GENOMES:
            sequences[name][column] = ancestral
        for isolate in node_isolates(node):
            sequences[isolate][column] = derived

    column, ancestral, derived, isolates = INCONGRUENT
    for name in GENOMES:
        sequences[name][column] = ancestral
    for isolate in isolates:
        sequences[isolate][column] = derived

    for name, base in zip(GENOMES, "AACGAA"):
        sequences[name][TRIALLELIC] = base
    for name, base in zip(GENOMES, ["A", "A", "A", "A", "A", "-"]):
        sequences[name][GAPPED] = base

    return {name: "".join(bases) for name, bases in sequences.items()}


def write_fasta(path: str, records: List[Tuple[str, str]], width: int = 70) -> None:
    with open(path, "w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), width):
                handle.write(sequence[start:start + width] + "\n")


@dataclass
class Dataset:
    """Paths to the written files, and what was planted in them."""

    directory: str
    alignment: str
    reference: str
    groups_csv: str
    groups_repgenr: str
    hc_yaml: str
    hc_repgenr: str
    gff: str
    trf: str
    sequences: Dict[str, str] = field(default_factory=dict)

    def planted(self) -> Dict[str, Tuple[int, str, str]]:
        return dict(PLANTED)

    def contig_of(self, column: int) -> Tuple[str, int]:
        """The contig and 1-based local position of an alignment column."""
        for name, start, end in CONTIGS:
            if start <= column < end:
                return name, column - start + 1
        raise AssertionError(f"column {column} is outside the reference")


@pytest.fixture
def dataset(tmp_path) -> Dataset:
    directory = str(tmp_path)
    sequences = build_alignment()

    alignment = os.path.join(directory, "alignment.fasta")
    write_fasta(alignment, [(name, sequences[name]) for name in GENOMES])

    reference = os.path.join(directory, "reference.fasta")
    write_fasta(
        reference,
        [(name, sequences["ref"][start:end]) for name, start, end in CONTIGS],
    )

    groups_csv = os.path.join(directory, "groups.csv")
    with open(groups_csv, "w") as handle:
        handle.write("Strain,Group\n")
        for cluster, isolates in CLUSTERS.items():
            for isolate in isolates:
                handle.write(f"{isolate},{cluster}\n")

    groups_repgenr = os.path.join(directory, "groups_repgenr.tsv")
    with open(groups_repgenr, "w") as handle:
        handle.write("representative\tcontained\tis_self\n")
        for cluster, isolates in CLUSTERS.items():
            handle.write(f"{cluster}\t{isolates[0]}.fasta\tTrue\n")
            for isolate in isolates[1:]:
                handle.write(f"{cluster}\t{isolate}.fasta\tFalse\n")

    hc_yaml = os.path.join(directory, "hc.yaml")
    with open(hc_yaml, "w") as handle:
        handle.write(
            "- root:\n"
            "  - A:\n"
            "    - A.1\n"
            "    - A.2\n"
            "  - B:\n"
            "    - B.1\n"
            "    - B.2\n"
        )

    hc_repgenr = os.path.join(directory, "hc_repgenr.tsv")
    with open(hc_repgenr, "w") as handle:
        handle.write("child\tparent\n")
        for child, parent in [
            ("A", "root"),
            ("B", "root"),
            ("A.1", "A"),
            ("A.2", "A"),
            ("B.1", "B"),
            ("B.2", "B"),
        ]:
            handle.write(f"{child}\t{parent}\n")

    gff = os.path.join(directory, "genes.gff")
    with open(gff, "w") as handle:
        handle.write("##gff-version 3\n")
        for name, start, end in CONTIGS:
            handle.write(
                f"{name}\tRefSeq\tregion\t1\t{end - start}\t.\t+\t.\tID={name}:1..{end - start}\n"
            )
            for contig, gene_start, gene_end, strand in GENES:
                if contig != name:
                    continue
                handle.write(
                    f"{name}\tRefSeq\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\tID=gene-{gene_start}\n"
                )

    trf = os.path.join(directory, "repeats.trf")
    with open(trf, "w") as handle:
        for name, _start, _end in CONTIGS:
            handle.write(f"@{name}\n")
            for contig, repeat_start, repeat_end in REPEATS:
                if contig == name:
                    handle.write(f"{repeat_start} {repeat_end} 4 8.0 4 90 0 50 25 25 25 25 1.9 ACGT ACGTACGT\n")

    dataset = Dataset(
        directory=directory,
        alignment=alignment,
        reference=reference,
        groups_csv=groups_csv,
        groups_repgenr=groups_repgenr,
        hc_yaml=hc_yaml,
        hc_repgenr=hc_repgenr,
        gff=gff,
        trf=trf,
        sequences=sequences,
    )
    return dataset


@pytest.fixture
def cansnps(dataset):
    """The extractor's output for the dataset: [position, node, base, base]."""
    from cansnpscorer.helpers.data_utils import canSNP_extractor, snp_extractor
    from cansnpscorer.helpers.file_utils import HCParser, fasta_parser, read_hc_file

    alignment = fasta_parser(dataset.alignment, "dict")
    snps = snp_extractor(alignment)
    tree = read_hc_file(dataset.hc_yaml, "yaml")
    groups = HCParser(dataset.groups_csv, data_format="csv").reader(tree)
    return [
        [int(position), node, *bases]
        for position, node, *bases in canSNP_extractor(alignment, groups, snps)
    ]


@pytest.fixture
def scored(dataset, cansnps):
    """The scored canSNPs of the dataset, as the evaluator computes them."""
    from cansnpscorer.cansnp_evaluator import annotations_by_contig, locate_snps
    from cansnpscorer.helpers.file_utils import gff_parser, reference_parser, trf_parser
    from cansnpscorer.helpers.scoring_methods import contexts_for, score_snp

    reference = reference_parser(dataset.reference)
    genes = annotations_by_contig(gff_parser(dataset.gff), reference, "GFF")
    repeats = annotations_by_contig(trf_parser(dataset.trf), reference, "TRF")
    contexts = contexts_for(locate_snps(cansnps, reference), genes, repeats)
    return [score_snp(context) for context in contexts if context.node != "None"]
