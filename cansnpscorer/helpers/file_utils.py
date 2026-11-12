
import csv
import logging
from collections import Counter
from dataclasses import dataclass
from operator import itemgetter
from typing import Dict, List, Optional, Sequence, Tuple, Union

from ruamel.yaml import YAML

from cansnpscorer.helpers.names import clean_name
from cansnpscorer.helpers.scoring_methods import TERMS, ScoredSnp
from cansnpscorer.hierarchy import Hierarchy, HierarchyError
from cansnpscorer.model import Reference


def read_hc_file(file_path: str, file_format: Optional[str] = None) -> Hierarchy:
    """
    Read the hierarchy of a scheme, whichever of the supported shapes it is in.

    Parameters:
        file_path: Path to the hierarchy file.
        file_format: One of 'yaml', 'json', 'newick' or 'repgenr'. Recognized
            from the file itself when not given.

    Returns:
        The hierarchy as a rooted tree.

    Raises:
        HierarchyError: If the format is not supported, or the file does not
            describe a single rooted tree.
        OSError: If the file cannot be read.
    """
    if file_format is None:
        file_format = recognize_hc_file_format(file_path)

    logging.info(f"Reading the {file_format} hierarchy in {file_path}")

    if file_format == 'yaml':
        yaml = YAML(typ='safe')
        with open(file_path) as handle:
            documents = list(yaml.load_all(handle))
        if not documents or documents[0] is None:
            raise HierarchyError(f"{file_path}: holds no YAML document")
        hierarchy = Hierarchy.from_nested(documents[0])

    elif file_format == 'json':
        with open(file_path) as handle:
            hierarchy = Hierarchy.from_json(handle.read())

    elif file_format == 'newick':
        with open(file_path) as handle:
            hierarchy = Hierarchy.from_newick(handle.read())

    elif file_format == 'repgenr':
        with open(file_path) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if not reader.fieldnames or not {'child', 'parent'} <= set(reader.fieldnames):
                raise HierarchyError(
                    f"{file_path}: a repgenr hierarchy needs 'child' and 'parent' columns, "
                    f"found {reader.fieldnames}")
            pairs = [(row['child'], row['parent']) for row in reader if row.get('child')]
        if not pairs:
            raise HierarchyError(f"{file_path}: holds no hierarchy rows")
        hierarchy = Hierarchy.from_child_parent(pairs)

    else:
        raise HierarchyError(
            f"{file_path}: unsupported hierarchy format '{file_format}'; "
            f"use YAML, JSON, Newick or the repgenr child and parent table")

    logging.info(
        f"Read a hierarchy of {len(hierarchy)} nodes rooted at '{hierarchy.root}', "
        f"with {len(hierarchy.leaves)} leaves")
    return hierarchy


def read_fasta(fasta_file: str) -> List[Tuple[str, str]]:
    """
    Read a FASTA file into (name, sequence) pairs, in file order.

    Sequence names are cleaned of non-alphanumeric characters, so that they
    match the strain names in a group file, and sequences are upper cased, so
    that soft masked regions are not mistaken for unusable bases.

    Args:
        fasta_file: Path to the FASTA file.

    Returns:
        The records in file order.

    Raises:
        ValueError: If the file holds no records, has sequence before its first
            header, or names a sequence more than once.
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
                name = clean_name(line[1:].strip())
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
            f"{fasta_file}: names {', '.join(repeated)} more than once; "
            "names must be unique once cleaned of non-alphanumeric characters")

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




class HCParser:
    """
    The isolates of each cluster, read from a group file.

    A group file names the leaves of the scheme and the isolates in each of
    them. How those leaves nest is the hierarchy's business, and combining the
    two is ``Hierarchy.isolates``.
    """

    def __init__(self, group_data: str, data_format: str = 'csv') -> None:
        """
        Read a group file.

        Args:
            group_data: The path to the file naming the clusters and their isolates.
            data_format: 'csv' or 'repgenr_group'.
        """
        logging.info(f"Reading the {data_format} group file {group_data}")
        self.group_label = "Group"
        self.strain_label = "Strain"
        self.isolates: Dict[str, List[str]] = {}
        if data_format == 'csv':
            self._group_handler(group_data)
        elif data_format == 'repgenr_group':
            self._repgenr_group_handler(group_data)
        else:
            raise ValueError(
                f"unknown group file format '{data_format}'; use 'csv' or 'repgenr_group'")
        logging.info(
            f"Read {sum(len(names) for names in self.isolates.values())} isolates in "
            f"{len(self.isolates)} clusters")

    def _repgenr_group_handler(self, repgenr_group_data_file: str) -> None:
        """
        Read a repgenr group file, in which each row places one isolate in the
        cluster of its representative.
        """
        with open(repgenr_group_data_file) as file:
            reader = csv.DictReader(file, delimiter='\t')
            required = {'representative', 'contained', 'is_self'}
            if not reader.fieldnames or not required <= set(reader.fieldnames):
                raise ValueError(
                    f"{repgenr_group_data_file}: a repgenr group file needs the columns "
                    f"{', '.join(sorted(required))}, found {reader.fieldnames}")
            for row in reader:
                # The representative names a cluster, which must read exactly as
                # the hierarchy names it; the contained value names a genome,
                # which must match the cleaned names of the alignment.
                representative = row['representative'].strip()
                contained = clean_name(row['contained'].replace('.fasta', ''))
                self.isolates.setdefault(representative, [])
                if contained not in self.isolates[representative]:
                    self.isolates[representative].append(contained)

    def _group_handler(self, csv_file: str) -> None:
        """
        Read a group file with Strain and Group columns, in either order.
        """
        with open(csv_file) as csvfile:
            reader = csv.reader(csvfile)
            for index, row in enumerate(reader):
                if not row:
                    continue
                if index == 0:
                    if self.group_label not in row or self.strain_label not in row:
                        raise ValueError(
                            f"{csv_file}: needs '{self.strain_label}' and '{self.group_label}' "
                            f"columns, found {row}")
                    self.group_index = row.index(self.group_label)
                    self.strain_index = row.index(self.strain_label)
                    continue
                if len(row) <= max(self.group_index, self.strain_index):
                    raise ValueError(
                        f"{csv_file}: row {index + 1} has {len(row)} fields, too few to hold "
                        f"both a strain and a group")
                group = row[self.group_index].strip()
                strain = clean_name(row[self.strain_index])
                if not group or not strain:
                    raise ValueError(
                        f"{csv_file}: row {index + 1} leaves the strain or the group empty")
                self.isolates.setdefault(group, []).append(strain)


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

#: The GFF feature types read as a reading frame. A gene and its CDS normally
#: coincide in bacteria, and 'gene' is what this tool has always read.
GENE_FEATURES = ("gene",)

#: GFF columns, by position, as the specification fixes them.
GFF_SEQID, GFF_TYPE, GFF_START, GFF_END, GFF_STRAND = 0, 2, 3, 4, 6


def gff_parser(gff_file: str) -> Dict[str, List[Tuple[int, int, str]]]:
    """
    Read the genes of a GFF file, keyed by the sequence they sit on.

    Columns are read by position, as the GFF specification fixes them, rather
    than by looking for a field whose value happens to be 'gene'. Only a line
    that begins with '#' is a comment, so a '#' inside an attribute no longer
    discards the feature. Sequences are keyed by their own names, cleaned the
    same way the FASTA readers clean them, so annotations are matched to
    contigs by name rather than by the order regions happen to appear in.

    Args:
        gff_file: The path to the GFF file.

    Returns:
        Gene intervals per sequence as (start, end, strand), 0-based and
        inclusive, sorted by start.

    Raises:
        OSError: If the file cannot be read.
    """
    logging.info(f"Reading genes from {gff_file}")
    genes: Dict[str, List[Tuple[int, int, str]]] = {}
    unreadable = 0

    with open(gff_file) as gff:
        for line_number, line in enumerate(gff, 1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= GFF_STRAND:
                continue
            if fields[GFF_TYPE] not in GENE_FEATURES:
                continue

            try:
                start = int(fields[GFF_START]) - 1
                end = int(fields[GFF_END]) - 1
            except ValueError:
                unreadable += 1
                logging.error(f"{gff_file}: line {line_number} has no readable start and end")
                continue

            if end < start:
                unreadable += 1
                logging.error(f"{gff_file}: line {line_number} ends before it starts")
                continue

            sequence = clean_name(fields[GFF_SEQID])
            strand = fields[GFF_STRAND] if fields[GFF_STRAND] in {"+", "-"} else "+"
            genes.setdefault(sequence, []).append((start, end, strand))

    for sequence in genes:
        genes[sequence].sort(key=itemgetter(0))

    if unreadable:
        logging.warning(f"Skipped {unreadable} unreadable features in {gff_file}")
    logging.info(
        f"Read {sum(len(found) for found in genes.values())} genes over "
        f"{len(genes)} sequences from {gff_file}")
    return genes


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
                sequence = clean_name(line[1:].split()[0]) if line[1:].split() else None
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

def recognize_group_file_format(file_path: str) -> str:
    """
    Recognize the file format by inspecting the first few lines.

    Args:
        file_path (str): The path to the file to inspect.

    Returns:
        str: The recognized file format ('csv' or 'repgenr_group').
    """
    with open(file_path) as file:
        lines = [line.strip() for line in file.readlines()[:5]]

    for line in lines:
        if {'Strain', 'Group'} <= set(field.strip() for field in line.split(',')):
            logging.info(f"{file_path} is a Strain and Group table")
            return 'csv'
        if {'representative', 'contained', 'is_self'} <= set(line.split('\t')):
            logging.info(f"{file_path} is a repgenr group file")
            return 'repgenr_group'

    raise ValueError(
        f"{file_path}: cannot tell which group file format this is; it needs either "
        f"'Strain' and 'Group' columns separated by commas, or the repgenr columns "
        f"'representative', 'contained' and 'is_self' separated by tabs")

def recognize_hc_file_format(file_path: str) -> str:
    """
    Recognize which of the supported hierarchy formats a file is in.

    The checks are ordered so that a format is only guessed once the more
    specific ones have been ruled out: a repgenr table is recognized by its
    columns, and a Newick string by its structure rather than by a stray
    semicolon.

    Returns:
        'repgenr', 'json', 'newick' or 'yaml'.

    Raises:
        HierarchyError: If the file is empty or matches none of them.
    """
    logging.info(f"Recognizing the format of {file_path}")

    with open(file_path) as handle:
        lines = [line.strip() for line in handle.readlines()[:20]]

    content = [line for line in lines if line and not line.startswith("#")]
    if not content:
        raise HierarchyError(f"{file_path}: is empty, or holds only comments")
    first_line = content[0]

    if "\t" in first_line and "child" in first_line and "parent" in first_line:
        recognized = "repgenr"
    elif first_line.startswith(("{", "[")):
        recognized = "json"
    elif ";" in "".join(content) and "(" in "".join(content):
        recognized = "newick"
    elif ":" in first_line or first_line.startswith("- "):
        recognized = "yaml"
    else:
        raise HierarchyError(
            f"{file_path}: cannot tell which hierarchy format this is from its first line "
            f"'{first_line}'; supported are YAML, JSON, Newick and the repgenr child and "
            f"parent table")

    logging.info(f"{file_path} is in {recognized} format")
    return recognized


