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

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

GENOMES = ["ref", "s1", "s1b", "s2", "s3", "s4"]
CLUSTERS = {"A.1": ["s1", "s1b"], "A.2": ["s2"], "B.1": ["s3"], "B.2": ["s4"]}
# The hierarchy as child and parent; the root is named only as a parent.
EDGES = [("A", "root"), ("B", "root"), ("A.1", "A"), ("A.2", "A"), ("B.1", "B"), ("B.2", "B")]

ALIGNMENT_LENGTH = 1000
CONTIGS = [("c1", 0, 700), ("c2", 700, 1000)]

# node -> (alignment column, ancestral base, derived base)
PLANTED: dict[str, tuple[int, str, str]] = {
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
# Tandem repeats: 1-based and inclusive, local to the contig.
REPEATS = [("c1", 60, 95), ("c2", 120, 160)]


def node_isolates(node: str) -> list[str]:
    """The isolates below a node, which is what a canSNP for it must separate."""
    descendants = {
        "root": ["A.1", "A.2", "B.1", "B.2"],
        "A": ["A.1", "A.2"],
        "B": ["B.1", "B.2"],
    }
    clusters = descendants.get(node, [node])
    return [isolate for cluster in clusters for isolate in CLUSTERS[cluster]]


def build_alignment() -> dict[str, str]:
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

    for name, base in zip(GENOMES, "AACGAA", strict=True):
        sequences[name][TRIALLELIC] = base
    for name, base in zip(GENOMES, ["A", "A", "A", "A", "A", "-"], strict=True):
        sequences[name][GAPPED] = base

    return {name: "".join(bases) for name, bases in sequences.items()}


def write_fasta(path: str, records: list[tuple[str, str]], width: int = 70) -> None:
    with open(path, "w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            handle.writelines(sequence[start:start + width] + "\n" for start in range(0, len(sequence), width))


@dataclass
class Dataset:
    """Paths to the written files, and what was planted in them."""

    directory: str
    alignment: str
    reference: str
    hierarchy: str
    hierarchy_newick: str
    genome_map: str
    gff: str
    repeats: str
    sequences: dict[str, str] = field(default_factory=dict)

    def planted(self) -> dict[str, tuple[int, str, str]]:
        return dict(PLANTED)

    def contig_of(self, column: int) -> tuple[str, int]:
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
    # Headers carry a description after the ID, as NCBI's do; the ID alone
    # names the contig.
    write_fasta(
        reference,
        [(f"{name} Testus organismus {name}, complete sequence", sequences["ref"][start:end])
         for name, start, end in CONTIGS],
    )

    # The genome map places each genome at a node, as FlexTaxD's
    # --genomeid2taxid and RepGenR2's genomes_map.tsv do: no header.
    genome_map = os.path.join(directory, "genome_map.tsv")
    with open(genome_map, "w") as handle:
        for cluster, isolates in CLUSTERS.items():
            handle.writelines(f"{isolate}\t{cluster}\n" for isolate in isolates)

    hierarchy = os.path.join(directory, "hierarchy.tsv")
    with open(hierarchy, "w") as handle:
        handle.write("child\tparent\n")
        handle.writelines(f"{child}\t{parent}\n" for child, parent in EDGES)

    hierarchy_newick = os.path.join(directory, "hierarchy.nwk")
    with open(hierarchy_newick, "w") as handle:
        handle.write("((A.1,A.2)A,(B.1,B.2)B)root;\n")

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
                    f"{name}\tRefSeq\tCDS\t{gene_start}\t{gene_end}\t.\t{strand}\t0\tID=cds-{gene_start};Parent=gene-{gene_start}\n"
                )

    # Tandem repeats as BED, 0-based and half open, as the trf subcommand
    # writes them.
    repeats = os.path.join(directory, "repeats.bed")
    with open(repeats, "w") as handle:
        for contig, repeat_start, repeat_end in REPEATS:
            handle.write(f"{contig}\t{repeat_start - 1}\t{repeat_end}\tperiod=4,copies=9.0\n")

    dataset = Dataset(
        directory=directory,
        alignment=alignment,
        reference=reference,
        hierarchy=hierarchy,
        hierarchy_newick=hierarchy_newick,
        genome_map=genome_map,
        gff=gff,
        repeats=repeats,
        sequences=sequences,
    )
    return dataset


@pytest.fixture
def isolates_by_node(dataset):
    """The isolates below each node of the dataset's hierarchy."""
    from cansnpscorer.helpers.file_utils import read_genome_map, read_hierarchy

    hierarchy = read_hierarchy(dataset.hierarchy)
    return hierarchy.isolates(read_genome_map(dataset.genome_map))


@pytest.fixture
def cansnps(dataset, isolates_by_node):
    """The canSNPs of the dataset, keyed by alignment column."""
    from cansnpscorer.helpers.data_utils import canSNP_extractor, snp_extractor
    from cansnpscorer.helpers.file_utils import read_alignment

    alignment = read_alignment(dataset.alignment)
    return canSNP_extractor(alignment, isolates_by_node, snp_extractor(alignment))


@pytest.fixture
def parsed_snps(dataset, cansnps):
    """The canSNPs as the SNP file records them, by contig and 1-based position."""
    from cansnpscorer.vcf import Snp

    parsed = []
    for cansnp in cansnps:
        contig, position = dataset.contig_of(cansnp.column)
        derived, ancestral = ((cansnp.derived, cansnp.ancestral) if cansnp.is_canonical
                              else cansnp.bases)
        parsed.append(Snp(contig, position, cansnp.node, derived, ancestral))
    return parsed


@pytest.fixture
def scored(dataset, parsed_snps):
    """The scored canSNPs of the dataset, as the evaluator computes them."""
    from cansnpscorer.cansnp_evaluator import annotations_by_contig, locate_snps
    from cansnpscorer.helpers.file_utils import read_gff3, read_intervals, reference_parser
    from cansnpscorer.helpers.scoring_methods import contexts_for, score_snp

    reference = reference_parser(dataset.reference)
    model = read_gff3(dataset.gff)
    genes = annotations_by_contig(model.genes, reference, "GFF")
    cds = annotations_by_contig(model.cds, reference, "CDS")
    repeats = annotations_by_contig(read_intervals(dataset.repeats), reference, "repeat")
    contexts = contexts_for(locate_snps(parsed_snps, reference), genes, repeats, cds)
    return [score_snp(context) for context in contexts if context.node != "None"]


# -- a scheme written by the export subcommand -------------------------------
#
# Shared by the audit tests, which read the scheme back in and change one thing
# at a time.

METADATA = {
    "organism": "testus_organismus",
    "genome": "TEST1",
    "strain": "testus",
    "citation": "Sjodin2026",
    "genbank_id": "GCA_000007125.1",
    "refseq_id": "GCF_000007125.1",
    "assembly_name": "ASM712v1",
}

TREE = "testus_organismus_tree.txt"
SNPS = "testus_organismus_snp.txt"
REFERENCES = "testus_organismus_references.txt"


def cli_run(*arguments):
    argv = sys.argv
    try:
        sys.argv = ["CanSNPCreator", *arguments]
        from cansnpscorer.cansnp_scorer import main as cli

        return cli()
    finally:
        sys.argv = argv


@pytest.fixture
def written(dataset, tmp_path):
    """A scheme written from the planted dataset, and its directory."""
    snp_file = str(tmp_path / "snps.vcf")
    scores = str(tmp_path / "scores.vcf")
    outdir = str(tmp_path / "scheme")
    assert cli_run("extract", "--genome-map", dataset.genome_map, "--hierarchy", dataset.hierarchy,
                   "--alignment", dataset.alignment, "--reference", dataset.reference,
                   "-o", snp_file) == 0
    assert cli_run("evaluate", "--snps", snp_file, "--reference", dataset.reference,
                   "--genes", dataset.gff, "--repeats", dataset.repeats,
                   "--hierarchy", dataset.hierarchy, "-o", scores,
                   "--cutoff", "0") == 0
    from cansnpscorer.cansnp_export import main as export_main

    assert export_main(scores, dataset.hierarchy, dataset.reference, outdir, **METADATA) == 0
    return outdir, scores


def edit(directory, name, replace, with_text):
    path = os.path.join(directory, name)
    with open(path) as handle:
        text = handle.read()
    assert replace in text, replace
    with open(path, "w") as handle:
        handle.write(text.replace(replace, with_text))


# -- VCF helpers -------------------------------------------------------------


def vcf_records(path):
    """The records of a VCF file this tool wrote."""
    from cansnpscorer import vcf

    return vcf.read(path).records


def vcf_lines(path):
    """The header and the record lines of a VCF file, as two lists."""
    with open(path) as handle:
        lines = handle.read().splitlines()
    return ([line for line in lines if line.startswith("#")],
            [line for line in lines if not line.startswith("#")])


def vcf_rewrite(path, out, keep=lambda fields: True, change=lambda fields: fields):
    """Copy a VCF file, dropping or changing record lines given as field lists."""
    header, records = vcf_lines(path)
    kept = []
    for line in records:
        fields = line.split("\t")
        if keep(fields):
            kept.append("\t".join(change(fields)))
    with open(out, "w") as handle:
        handle.write("\n".join(header + kept) + "\n")
    return str(out)


def info_of(fields):
    """The INFO column of a record line as a dictionary."""
    return dict(item.partition("=")[::2] for item in fields[7].split(";"))
