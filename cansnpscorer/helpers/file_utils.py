
import csv
import logging
from collections import Counter
from dataclasses import dataclass
from operator import itemgetter
from typing import Dict, List, Optional, Sequence, Tuple, Union

from cansnpscorer.helpers.names import strip_fasta_suffix
from cansnpscorer.helpers.scoring_methods import TERMS, ScoredSnp
from cansnpscorer.hierarchy import Hierarchy, HierarchyError
from cansnpscorer.model import Reference

HIERARCHY_FORMATS = ("table", "newick")


def hierarchy_format(path: str) -> str:
    """Tell a Newick tree from a child and parent table by their first character.

    Returns:
        'newick' or 'table'.

    Raises:
        HierarchyError: If the file is empty or is neither.
    """
    with open(path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if stripped.startswith("("):
                return "newick"
            columns = stripped.split("\t")
            if "child" in columns and "parent" in columns:
                return "table"
            raise HierarchyError(
                f"{path}: is neither a Newick tree nor a tab separated table with 'child' and "
                f"'parent' columns; its first line is '{stripped[:60]}'")
    raise HierarchyError(f"{path}: is empty")


def read_hierarchy(path: str) -> Hierarchy:
    """
    Read the hierarchy of a scheme: which node lies under which.

    Two formats are read. A tab separated table with a header naming 'child'
    and 'parent' columns holds one row per node; the root's parent is empty,
    or the root is the one node named only as a parent. This is the layout of
    FlexTaxD's taxonomy and RepGenR2's tree2tax.tsv, and other columns, such as
    tree2tax's rank, are ignored. A Newick tree names its nodes in its labels;
    every label is read as a name, so a node called 2.1 stays 2.1.

    Raises:
        HierarchyError: If the file is in neither format, or the nodes do not
            form a single rooted tree.
        OSError: If the file cannot be read.
    """
    file_format = hierarchy_format(path)
    logging.info(f"Reading the hierarchy in {path} as a {file_format}")

    if file_format == "newick":
        with open(path) as handle:
            hierarchy = Hierarchy.from_newick(handle.read())
    else:
        pairs = []
        with open(path, newline="") as handle:
            lines = [line.rstrip("\r\n") for line in handle
                     if line.strip() and not line.startswith("#")]
        header = lines[0].split("\t")
        child_at, parent_at = header.index("child"), header.index("parent")
        for number, line in enumerate(lines[1:], 2):
            fields = line.split("\t")
            if len(fields) <= child_at or not fields[child_at].strip():
                raise HierarchyError(f"{path}: row {number} names no child")
            parent = fields[parent_at] if len(fields) > parent_at else ""
            pairs.append((fields[child_at], parent))
        if not pairs:
            raise HierarchyError(f"{path}: holds no rows below its header")
        hierarchy = Hierarchy.from_child_parent(pairs)

    logging.info(
        f"Read a hierarchy of {len(hierarchy)} nodes rooted at '{hierarchy.root}', "
        f"with {len(hierarchy.leaves)} leaves")
    return hierarchy


def read_fasta(fasta_file: str) -> List[Tuple[str, str]]:
    """
    Read a FASTA file into (name, sequence) pairs, in file order.

    A record is named by its ID, the first word of its header, as GFF3, VCF,
    BED and samtools faidx all name a sequence; the rest of the header is a
    description. Sequences are upper cased, so that soft masked regions are
    not mistaken for unusable bases.

    Args:
        fasta_file: Path to the FASTA file.

    Returns:
        The records in file order.

    Raises:
        ValueError: If the file holds no records, has sequence before its first
            header, has a header without an ID, or names a sequence more than once.
    """
    logging.info(f"Reading FASTA records from {fasta_file}")

    records: List[Tuple[str, str]] = []
    name: Optional[str] = None
    chunks: List[str] = []
    masked = 0

    with open(fasta_file) as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                words = line[1:].split()
                if not words:
                    raise ValueError(f"{fasta_file}: line {line_number} is a header without an ID")
                name = words[0]
                chunks = []
            else:
                if name is None:
                    raise ValueError(
                        f"{fasta_file}: sequence before the first header, at line {line_number}")
                masked += sum(1 for base in line if base.islower())
                chunks.append(line.upper())

    if name is not None:
        records.append((name, "".join(chunks)))

    if not records:
        raise ValueError(f"{fasta_file}: holds no FASTA records")

    counts = Counter(record_name for record_name, _ in records)
    repeated = sorted(record_name for record_name, count in counts.items() if count > 1)
    if repeated:
        raise ValueError(
            f"{fasta_file}: names {', '.join(repeated)} more than once; sequence IDs, the first "
            "word of each header, must be unique")

    if masked:
        logging.info(f"Upper cased {masked} soft masked bases from {fasta_file}")
    logging.info(f"Read {len(records)} records from {fasta_file}")
    return records


def fasta_parser(fasta_file: str, output_format: str = 'list') -> Union[List[str], List[Tuple[str, str]], Dict[str, str]]:
    """
    Parse a FASTA file.

    Args:
        fasta_file: Path to the FASTA file.
        output_format: Output format ('list', 'dict', 'named_list').

    Returns:
        Parsed genomes in the requested output format.
    """
    records = read_fasta(fasta_file)

    if output_format == 'dict':
        return dict(records)
    if output_format == 'named_list':
        return records
    return [sequence for _name, sequence in records]


def read_alignment(alignment_file: str) -> Dict[str, str]:
    """
    Read a whole genome alignment, with the reference genome as its first record.

    Args:
        alignment_file: Path to the alignment in FASTA format.

    Returns:
        The aligned sequences, keyed by name, with the reference first.

    Raises:
        ValueError: If the alignment holds fewer than two records, or if the
            records are not all of the same length.
    """
    records = read_fasta(alignment_file)

    if len(records) < 2:
        raise ValueError(
            f"{alignment_file}: an alignment needs the reference and at least one other "
            f"genome, found {len(records)}")

    lengths = {len(sequence) for _name, sequence in records}
    if len(lengths) > 1:
        shortest = min(records, key=lambda record: len(record[1]))
        longest = max(records, key=lambda record: len(record[1]))
        raise ValueError(
            f"{alignment_file}: the records are not aligned; '{longest[0]}' holds "
            f"{len(longest[1])} columns and '{shortest[0]}' holds {len(shortest[1])}")

    logging.info(
        f"Read an alignment of {len(records)} genomes over {len(records[0][1])} columns, "
        f"with '{records[0][0]}' as the reference")
    return dict(records)




def read_genome_map(path: str) -> Dict[str, List[str]]:
    """
    Read which genomes sit at each node: two tab separated columns, genome and
    node, and no header.

    This is the layout of FlexTaxD's --genomeid2taxid file and of RepGenR2's
    genomes_map.tsv. A genome may sit at an internal node. Genome names are
    kept as written and matched to alignment records later, so a name may be
    a file name, a record ID or an assembly accession. Lines beginning with
    '#' are comments.

    Returns:
        The genomes at each node, in file order.

    Raises:
        ValueError: If a line does not hold exactly a genome and a node, or a
            genome is placed at two nodes.
        OSError: If the file cannot be read.
    """
    logging.info(f"Reading the genome map {path}")
    clusters: Dict[str, List[str]] = {}
    placed: Dict[str, str] = {}
    with open(path, newline="") as handle:
        for number, line in enumerate(handle, 1):
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 2 or not fields[0].strip() or not fields[1].strip():
                raise ValueError(
                    f"{path}: line {number} should hold a genome and a node separated by one "
                    f"tab, with no header line")
            genome = strip_fasta_suffix(fields[0].strip())
            node = fields[1].strip()
            if genome in placed and placed[genome] != node:
                raise ValueError(
                    f"{path}: line {number} places {genome} at {node}, but it is already "
                    f"placed at {placed[genome]}")
            if genome not in placed:
                placed[genome] = node
                clusters.setdefault(node, []).append(genome)
    if not clusters:
        raise ValueError(f"{path}: places no genome")
    logging.info(f"Read {len(placed)} genomes placed at {len(clusters)} nodes")
    return clusters


#: The columns of the SNP file the extractor writes.
SNP_HEADER = ("contig", "position", "node", "derived_base", "ancestral_base")


@dataclass
class ParsedSnp:
    """One row of a SNP file.

    A SNP from the current format is addressed by contig and 1-based position.
    One from the older format carries a position in the concatenated genome
    instead, counted from zero, and no polarity.
    """

    node: str
    alleles: Tuple[str, ...]
    contig: Optional[str] = None
    position: Optional[int] = None
    global_position: Optional[int] = None

    @property
    def derived(self) -> Optional[str]:
        return self.alleles[0] if self.contig is not None and self.alleles else None

    @property
    def ancestral(self) -> Optional[str]:
        return self.alleles[1] if self.contig is not None and len(self.alleles) > 1 else None


def snp_file_parser(snpfile: str) -> List[ParsedSnp]:
    """
    Read a SNP file written by the extract subcommand.

    Both the current tab separated format, which names the contig and gives a
    1-based position, and the older space separated format, which gives a
    position in the concatenated genome counted from zero, are accepted. The
    older format is reported once as what it is, since its positions cannot be
    checked against a contig and carry no polarity.

    Args:
        snpfile: The path to the SNP file.

    Returns:
        The SNPs in file order.

    Raises:
        ValueError: If the file holds no readable rows.
        OSError: If the file cannot be read.
    """
    logging.info(f"Reading SNPs from {snpfile}")

    snps: List[ParsedSnp] = []
    legacy = 0
    unreadable = 0

    with open(snpfile) as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            if line.split("\t")[0:1] == [SNP_HEADER[0]]:
                continue

            fields = line.split("\t") if "\t" in line else line.split(" ")
            fields = [field for field in fields if field]

            try:
                if "\t" in line:
                    contig, position, node, *bases = fields
                    snps.append(ParsedSnp(node=node, alleles=tuple(bases), contig=contig,
                                          position=int(position)))
                else:
                    position, node, *bases = fields
                    legacy += 1
                    snps.append(ParsedSnp(node=node, alleles=tuple(bases),
                                          global_position=int(position)))
            except ValueError:
                unreadable += 1
                logging.error(f"{snpfile}: cannot read line {line_number}: {line}")

    if unreadable:
        logging.warning(f"Skipped {unreadable} unreadable rows in {snpfile}")
    if legacy:
        logging.warning(
            f"{legacy} rows of {snpfile} are in the older format, which counts positions over "
            "the whole genome from zero and records no derived or ancestral base. Re-run "
            "extract to get positions by contig.")
    if not snps:
        raise ValueError(f"{snpfile}: holds no readable SNP rows")

    logging.info(f"Read {len(snps)} SNPs from {snpfile}")
    return snps

def reference_parser(reference: str) -> Reference:
    """
    Parse a reference genome file into contigs.

    Args:
        reference: The path to the reference genome in FASTA format.

    Returns:
        The reference, with one contig per record in file order.

    Raises:
        ValueError: If the file is not usable as a reference.
    """
    logging.info(f"Reading the reference genome from {reference}")
    parsed = Reference.from_records(read_fasta(reference))
    logging.info(
        f"Read {len(parsed)} contigs totalling {parsed.length} bases: {', '.join(parsed.names)}")
    return parsed

#: GFF3 feature types that place a SNP inside a gene. A CDS counts as well as a
#: gene, since Prokka writes CDS features without gene features by default.
GENE_FEATURES = frozenset({"gene", "CDS"})

GFF_SEQID, GFF_TYPE, GFF_START, GFF_END, GFF_STRAND, GFF_PHASE = 0, 2, 3, 4, 6, 7


@dataclass
class GeneModel:
    """The genes of a GFF3 file, by sequence ID, 0-based and inclusive.

    ``genes`` holds every gene and CDS as (start, end, strand) and decides
    whether a SNP lies in a gene. ``cds`` holds the CDS features as (start,
    end, strand, phase), which give the reading frame a codon is read in.
    """

    genes: Dict[str, List[Tuple[int, int, str]]]
    cds: Dict[str, List[Tuple[int, int, str, int]]]


def read_gff3(gff_file: str) -> GeneModel:
    """
    Read the genes and coding sequences of a GFF3 file.

    Columns are read by position, as GFF3 fixes them, and sequences are named
    by their seqid exactly, as the reference FASTA names them by its IDs. A
    line beginning with '#' is a comment or directive, and a ##FASTA directive
    ends the annotation.

    Raises:
        ValueError: If a feature has no readable start and end.
        OSError: If the file cannot be read.
    """
    logging.info(f"Reading genes from {gff_file}")
    genes: Dict[str, List[Tuple[int, int, str]]] = {}
    cds: Dict[str, List[Tuple[int, int, str, int]]] = {}

    with open(gff_file) as gff:
        for line_number, line in enumerate(gff, 1):
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"{gff_file}: line {line_number} holds {len(fields)} tab separated columns "
                    f"where GFF3 has 9")
            if fields[GFF_TYPE] not in GENE_FEATURES:
                continue
            try:
                start = int(fields[GFF_START]) - 1
                end = int(fields[GFF_END]) - 1
            except ValueError:
                raise ValueError(
                    f"{gff_file}: line {line_number} has no readable start and end") from None
            if end < start or start < 0:
                raise ValueError(f"{gff_file}: line {line_number} has an impossible interval")
            sequence = fields[GFF_SEQID]
            strand = fields[GFF_STRAND] if fields[GFF_STRAND] in {"+", "-"} else "+"
            genes.setdefault(sequence, []).append((start, end, strand))
            if fields[GFF_TYPE] == "CDS":
                phase = fields[GFF_PHASE]
                if phase not in {"0", "1", "2"}:
                    raise ValueError(
                        f"{gff_file}: line {line_number} is a CDS without a phase of 0, 1 or 2")
                cds.setdefault(sequence, []).append((start, end, strand, int(phase)))

    # A gene and its CDS usually cover the same bases; one interval is enough
    # to say a SNP is in a gene.
    for sequence in genes:
        genes[sequence] = sorted(set(genes[sequence]), key=itemgetter(0, 1))
    for sequence in cds:
        cds[sequence].sort(key=itemgetter(0, 1))
    logging.info(
        f"Read {sum(len(found) for found in genes.values())} genes and CDS, "
        f"{sum(len(found) for found in cds.values())} of them CDS, over {len(genes)} "
        f"sequences from {gff_file}")
    return GeneModel(genes=genes, cds=cds)


def trf_parser(trf_file_path: str) -> Dict[str, List[Tuple[int, int]]]:
    """
    Read the tandem repeats of a Tandem Repeats Finder -ngs file.

    Args:
        trf_file_path: The path to the TRF output.

    Returns:
        Repeat intervals per sequence as (start, end), 0-based and inclusive,
        sorted by start. TRF reports them 1-based, as a GFF does.

    Raises:
        OSError: If the file cannot be read.
    """
    logging.info(f"Reading tandem repeats from {trf_file_path}")
    intervals: Dict[str, List[Tuple[int, int]]] = {}
    sequence: Optional[str] = None
    orphaned = 0

    with open(trf_file_path) as trf_file:
        for line in trf_file:
            line = line.strip()
            if not line:
                continue
            if line.startswith("@"):
                sequence = line[1:].split()[0] if line[1:].split() else None
                if sequence:
                    intervals.setdefault(sequence, [])
                continue

            fields = line.split()
            try:
                start, end = int(fields[0]) - 1, int(fields[1]) - 1
            except (IndexError, ValueError):
                continue
            if sequence is None:
                orphaned += 1
                continue
            intervals[sequence].append((start, end) if start <= end else (end, start))

    for name in intervals:
        intervals[name].sort(key=itemgetter(0))

    if orphaned:
        logging.warning(
            f"{trf_file_path}: {orphaned} repeats appear before any sequence header and were "
            "left out")
    logging.info(
        f"Read {sum(len(found) for found in intervals.values())} repeats over "
        f"{len(intervals)} sequences from {trf_file_path}")
    return intervals


def score_file_writer(scores: Sequence[ScoredSnp], output_file_path: str) -> None:
    """
    Write scored canSNPs to a CSV file.

    Args:
        scores: The scored canSNPs to write.
        output_file_path: Where to write the file.
    """
    logging.info(f"Writing {len(scores)} scored canSNPs to {output_file_path}")

    with open(output_file_path, 'w', newline='') as out_file:
        csv_writer = csv.writer(out_file, delimiter=',')
        csv_writer.writerow(ScoredSnp.HEADER)
        csv_writer.writerows(score.row() for score in scores)


def read_score_file(path: str) -> List[ScoredSnp]:
    """
    Read a score file written by the evaluate subcommand.

    Args:
        path: The path to the CSV file.

    Returns:
        The scored canSNPs in file order.

    Raises:
        ValueError: If the file does not carry the expected columns, or holds no rows.
        OSError: If the file cannot be read.
    """
    logging.info(f"Reading scores from {path}")
    scores: List[ScoredSnp] = []

    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        missing = [column for column in ScoredSnp.HEADER if column not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(
                f"{path}: is missing the columns {', '.join(missing)}; it should be the output "
                f"of the evaluate subcommand")

        for line_number, row in enumerate(reader, 2):
            try:
                scores.append(ScoredSnp(
                    contig=row["Contig"],
                    position=int(row["Position"]),
                    node=row["Group"],
                    total=int(row["Score"]),
                    terms={term: int(row[column]) for term, column
                           in zip(TERMS, ScoredSnp.HEADER[6:])},
                    derived=row["Derived"] or None,
                    ancestral=row["Ancestral"] or None,
                ))
            except (TypeError, ValueError):
                logging.error(f"{path}: cannot read line {line_number}")

    if not scores:
        raise ValueError(f"{path}: holds no readable score rows")

    logging.info(f"Read {len(scores)} scores from {path}")
    return scores


def cansnp_file_writer(rows: Sequence[Sequence[Union[str, int]]], output_filepath: str) -> None:
    """
    Write the extracted SNPs as a tab separated file with a header.

    Args:
        rows: One row per SNP, following SNP_HEADER.
        output_filepath: Where to write the file.
    """
    logging.info(f"Writing {len(rows)} SNPs to {output_filepath}")

    with open(output_filepath, 'w') as out_file:
        out_file.write("\t".join(SNP_HEADER) + "\n")
        for row in rows:
            out_file.write("\t".join("" if field is None else str(field) for field in row) + "\n")
