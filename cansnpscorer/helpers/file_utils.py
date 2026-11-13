
import gzip
import logging
from collections import Counter
from dataclasses import dataclass, field
from operator import itemgetter
from urllib.parse import unquote

from cansnpscorer.helpers.names import strip_fasta_suffix
from cansnpscorer.hierarchy import Hierarchy, HierarchyError
from cansnpscorer.model import Reference

HIERARCHY_FORMATS = ("table", "newick")


def open_text(path: str, newline: str | None = None):
    """Open a text file for reading, gzip or bgzip compressed or not.

    NCBI distributes genomes and annotations gzipped, and VCF is often
    bgzipped; both are recognised by their first two bytes.
    """
    with open(path, "rb") as probe:
        compressed = probe.read(2) == b"\x1f\x8b"
    # utf-8-sig drops a byte order mark, which would otherwise become part of
    # the first header or column name.
    if compressed:
        return gzip.open(path, "rt", newline=newline, encoding="utf-8-sig")
    return open(path, newline=newline, encoding="utf-8-sig")


def hierarchy_format(path: str) -> str:
    """Tell a Newick tree from a child and parent table by their first character.

    Returns:
        'newick' or 'table'.

    Raises:
        HierarchyError: If the file is empty or is neither.
    """
    with open_text(path) as handle:
        in_comment = False
        for line in handle:
            stripped = line.strip()
            # A Newick tree may open with comments in square brackets, as
            # [&R] marks a rooted tree; they are skipped.
            while in_comment or stripped.startswith("["):
                end = stripped.find("]")
                in_comment = end < 0
                stripped = "" if in_comment else stripped[end + 1:].lstrip()
                if in_comment:
                    break
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
    tree2tax's rank, are ignored. A root given as its own parent, as NCBI's
    nodes.dmp gives it, is read as the root. A Newick tree names its nodes in its labels;
    every label is read as a name, so a node called 2.1 stays 2.1.

    Raises:
        HierarchyError: If the file is in neither format, or the nodes do not
            form a single rooted tree.
        OSError: If the file cannot be read.
    """
    file_format = hierarchy_format(path)
    logging.info(f"Reading the hierarchy in {path} as a {file_format}")

    if file_format == "newick":
        with open_text(path) as handle:
            hierarchy = Hierarchy.from_newick(handle.read())
    else:
        pairs = []
        with open_text(path, newline="") as handle:
            lines = [line.rstrip("\r\n") for line in handle
                     if line.strip() and not line.startswith("#")]
        header = lines[0].split("\t")
        child_at, parent_at = header.index("child"), header.index("parent")
        for number, line in enumerate(lines[1:], 2):
            fields = line.split("\t")
            if len(fields) <= child_at or not fields[child_at].strip():
                raise HierarchyError(f"{path}: row {number} names no child")
            parent = fields[parent_at] if len(fields) > parent_at else ""
            if parent.strip() == fields[child_at].strip():
                parent = ""  # NCBI's nodes.dmp names the root as its own parent
            pairs.append((fields[child_at], parent))
        if not pairs:
            raise HierarchyError(f"{path}: holds no rows below its header")
        hierarchy = Hierarchy.from_child_parent(pairs)

    logging.info(
        f"Read a hierarchy of {len(hierarchy)} nodes rooted at '{hierarchy.root}', "
        f"with {len(hierarchy.leaves)} leaves")
    return hierarchy


#: IUPAC nucleotide codes, and the gap and unknown characters aligners write.
SEQUENCE_CHARACTERS = frozenset("ACGTURYSWKMBDHVN-.~*?")


def read_fasta(fasta_file: str) -> list[tuple[str, str]]:
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

    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    masked = 0

    with open_text(fasta_file) as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith(";"):
                continue  # ';' opens a comment line in the original FASTA format
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                words = line[1:].split()
                if not words or line[1:2].isspace():
                    raise ValueError(f"{fasta_file}: line {line_number} is a header without an ID "
                                     f"directly after '>'")
                name = words[0]
                chunks = []
            else:
                if name is None:
                    raise ValueError(
                        f"{fasta_file}: sequence before the first header, at line {line_number}")
                line = "".join(line.split())
                invalid = set(line.upper()) - SEQUENCE_CHARACTERS
                if invalid:
                    raise ValueError(
                        f"{fasta_file}: line {line_number} holds characters that are not IUPAC "
                        f"bases or gaps: {''.join(sorted(invalid))[:10]}")
                masked += sum(1 for base in line if base.islower())
                chunks.append(line.upper())

    if name is not None:
        records.append((name, "".join(chunks)))

    if not records:
        raise ValueError(f"{fasta_file}: holds no FASTA records")
    empty = [record_name for record_name, sequence in records if not sequence]
    if empty:
        raise ValueError(f"{fasta_file}: records hold no sequence: {', '.join(empty[:5])}")

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


def fasta_ids(path: str) -> list[str]:
    """The record IDs of a FASTA file, in order, without keeping the sequences."""
    ids: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                words = line[1:].split()
                if words:
                    ids.append(words[0])
    return ids


def read_alignment(alignment_file: str) -> dict[str, str]:
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




def read_genome_map(path: str) -> dict[str, list[str]]:
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
    clusters: dict[str, list[str]] = {}
    placed: dict[str, str] = {}
    with open_text(path, newline="") as handle:
        for number, line in enumerate(handle, 1):
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            if not placed and [field.strip().lower() for field in fields] == ["genome", "node"]:
                continue  # a header naming the columns, which the format does not need
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

#: Genetic codes whose codons translate as the standard table does. Table 11,
#: for bacteria, differs from table 1 only in its start codons.
STANDARD_CODES = frozenset({"1", "11"})

GFF_SEQID, GFF_TYPE, GFF_START, GFF_END, GFF_STRAND, GFF_PHASE, GFF_ATTRIBUTES = (
    0, 2, 3, 4, 6, 7, 8)


@dataclass
class GeneModel:
    """The genes of a GFF3 file, by sequence ID, 0-based and inclusive.

    ``genes`` holds the bases covered by any gene or CDS as (start, end),
    joined where they overlap, and decides whether a SNP lies in a gene.
    ``cds`` holds the CDS features that give a reading frame as (start, end,
    strand, phase). ``sequences`` holds the sequences a ##FASTA section
    carries, as Prokka and Bakta write them, by which the annotated sequences
    can be matched to the reference's contigs when their names differ.
    ``loci`` holds each gene on its own as (start, end), sorted and not
    joined, which tells one gene from the next where they abut or overlap.
    """

    genes: dict[str, list[tuple[int, int]]]
    cds: dict[str, list[tuple[int, int, str, int]]]
    sequences: dict[str, str] = field(default_factory=dict)
    loci: dict[str, list[tuple[int, int]]] = field(default_factory=dict)


def gff3_attributes(text: str) -> dict[str, str]:
    """The attributes of a GFF3 feature, percent-decoding their values."""
    attributes = {}
    for item in text.strip().split(";"):
        key, _equals, value = item.partition("=")
        if key:
            attributes[key.strip()] = unquote(value)
    return attributes


def read_gff3(gff_file: str) -> GeneModel:
    """
    Read the genes and coding sequences of a GFF3 file.

    Prokka's, Bakta's and NCBI's GFF3 are read alike. Columns are read by
    position, as GFF3 fixes them, and sequences are named by their seqid. A
    line beginning with '#' is a comment or directive; a ##FASTA directive ends
    the features, and the sequences after it are kept.

    A CDS gives no reading frame when it is a pseudogene's (NCBI marks it
    pseudo=true) or uses a genetic code whose codons translate differently from
    the standard code (transl_table other than 1 or 11). It still places the
    SNP in a gene.

    Raises:
        ValueError: If a feature is malformed.
        OSError: If the file cannot be read.
    """
    logging.info(f"Reading genes from {gff_file}")
    genes: dict[str, list[tuple[int, int, str]]] = {}
    features: dict[str, dict[str, set[tuple[int, int]]]] = {}
    cds: dict[str, list[tuple[int, int, str, int]]] = {}
    sequences: dict[str, list[str]] = {}
    pseudo = 0
    other_code: Counter = Counter()

    with open_text(gff_file) as gff:
        in_fasta = False
        name = None
        for line_number, line in enumerate(gff, 1):
            line = line.rstrip("\r\n")
            if in_fasta:
                if line.startswith(">"):
                    words = line[1:].split()
                    name = words[0] if words else None
                    if name is not None:
                        sequences[name] = []
                elif name is not None and line.strip():
                    sequences[name].append(line.strip().upper())
                continue
            if line.startswith("##FASTA"):
                in_fasta = True
                continue
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t")
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
            sequence = unquote(fields[GFF_SEQID])
            strand = fields[GFF_STRAND] if fields[GFF_STRAND] in {"+", "-"} else "+"
            genes.setdefault(sequence, []).append((start, end, strand))
            features.setdefault(sequence, {}).setdefault(fields[GFF_TYPE], set()).add(
                (start, end))
            if fields[GFF_TYPE] != "CDS":
                continue
            phase = fields[GFF_PHASE]
            if phase not in {"0", "1", "2"}:
                raise ValueError(
                    f"{gff_file}: line {line_number} is a CDS without a phase of 0, 1 or 2")
            attributes = gff3_attributes(fields[GFF_ATTRIBUTES])
            if attributes.get("pseudo", "").lower() == "true" or "pseudogene" in attributes:
                pseudo += 1
                continue
            code = attributes.get("transl_table", "11")
            if code not in STANDARD_CODES:
                other_code[code] += 1
                continue
            cds.setdefault(sequence, []).append((start, end, strand, int(phase)))

    # Whether a SNP is in a gene does not depend on the strand, and the genes
    # are searched with a binary search, which needs intervals that do not
    # overlap: a gene with its CDS, or a gene holding smaller features, would
    # otherwise hide the bases it covers.
    merged = {sequence: merge_intervals([(start, end) for start, end, _ in found])
              for sequence, found in genes.items()}
    for sequence in cds:
        cds[sequence].sort(key=itemgetter(0, 1))
    if pseudo:
        logging.info(f"{pseudo} pseudogene CDS of {gff_file} give no reading frame")
    for code, count in sorted(other_code.items()):
        logging.warning(
            f"{count} CDS of {gff_file} use genetic code {code}, whose codons translate "
            f"differently from the standard code, and give no reading frame")
    logging.info(
        f"Read {sum(len(found) for found in genes.values())} genes and CDS features, "
        f"{sum(len(found) for found in cds.values())} CDS with a reading frame, over "
        f"{len(genes)} sequences from {gff_file}"
        + (f", and {len(sequences)} sequences after ##FASTA" if sequences else ""))
    # Each gene on its own, for telling which gene a SNP lies in: the gene
    # features where the file has them, as NCBI and Bakta write a gene with
    # its CDS, and the CDS features where it does not, as Prokka can.
    loci = {sequence: sorted(found.get("gene") or found.get("CDS", ()))
            for sequence, found in features.items()}
    return GeneModel(genes=merged, cds=cds, loci=loci,
                     sequences={key: "".join(value) for key, value in sequences.items()})


def read_contig_aliases(path: str) -> list[set[str]]:
    """
    Read which names denote the same sequence.

    Two standard layouts are read. NCBI's assembly_report.txt, which sits
    beside every assembly NCBI distributes, gives each molecule's sequence
    name, GenBank accession, RefSeq accession and UCSC name; a RefSeq GFF3
    names the chromosome of Coxiella burnetii RSA 493 NC_002971.4 where the
    GenBank FASTA names it AE016828.3. A UCSC chromAlias table gives every
    name of one sequence on one tab separated line.

    Returns:
        One set of names per sequence.

    Raises:
        ValueError: If the file holds no names.
    """
    groups: list[set[str]] = []
    with open_text(path, newline="") as handle:
        lines = [line.rstrip("\r\n") for line in handle]
    ncbi = any(line.startswith("# Sequence-Name") for line in lines)
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if ncbi:
            if len(fields) < 7:
                continue
            names = [fields[index] for index in (0, 4, 6, 9) if index < len(fields)]
        else:
            names = fields
        group = {name.strip() for name in names if name.strip() and name.strip() != "na"}
        if len(group) > 1:
            groups.append(group)
    if not groups:
        raise ValueError(f"{path}: names no sequence by more than one name")
    logging.info(f"Read {len(groups)} sequences with more than one name from {path}")
    return groups


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Sort intervals and join those that overlap or touch."""
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def read_intervals(path: str) -> dict[str, list[tuple[int, int]]]:
    """
    Read regions of the reference from BED or GFF3: tandem repeats, or regions
    to mask such as the recombinant regions Gubbins reports.

    A file whose leading comments hold a ##gff-version directive is GFF3, and
    every feature in it is a region, 1-based and inclusive. Anything else is
    BED: chrom, start and end in the first three columns, 0-based and half
    open, with track, browser and '#' lines skipped. Sequences are named
    exactly as the reference FASTA's IDs name them.

    Returns:
        The regions per sequence as (start, end), 0-based and inclusive,
        sorted, with overlapping regions joined.

    Raises:
        ValueError: If a line is not a well formed region.
        OSError: If the file cannot be read.
    """
    regions: dict[str, list[tuple[int, int]]] = {}
    with open_text(path, newline="") as handle:
        lines = [line.rstrip("\r\n") for line in handle]
    is_gff = False
    for line in lines:
        if not line.strip():
            continue
        if line.startswith("##gff-version"):
            is_gff = True
            break
        if not line.startswith("#"):
            break

    for number, line in enumerate(lines, 1):
        if not line.strip():
            continue
        if is_gff:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"{path}: line {number} holds {len(fields)} tab separated columns where "
                    f"GFF3 has 9")
            sequence, first, last = fields[0], fields[3], fields[4]
            try:
                start, end = int(first) - 1, int(last) - 1
            except ValueError:
                raise ValueError(f"{path}: line {number} has no readable start and end") from None
        else:
            if line.startswith("#") or line.split()[0] in ("track", "browser"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                # BED is tab separated, but whitespace separated files exist.
                fields = line.split()
            if len(fields) < 3:
                raise ValueError(
                    f"{path}: line {number} holds {len(fields)} tab separated columns where BED "
                    f"has at least 3")
            sequence = fields[0]
            try:
                start, end = int(fields[1]), int(fields[2]) - 1
            except ValueError:
                raise ValueError(f"{path}: line {number} has no readable start and end") from None
            if end == start - 1 and start >= 0:
                continue  # BED allows an empty region, which covers no base
        if start < 0 or end < start:
            raise ValueError(f"{path}: line {number} has an impossible interval")
        regions.setdefault(sequence, []).append((start, end))

    merged = {sequence: merge_intervals(found) for sequence, found in regions.items()}
    logging.info(
        f"Read {sum(len(found) for found in merged.values())} regions over {len(merged)} "
        f"sequences from {path}, as {'GFF3' if is_gff else 'BED'}")
    return merged
