"""SNPs and their scores as VCF 4.3.

The extract subcommand writes every SNP of the alignment, and the evaluate
subcommand writes the canonical ones with their scores, both as VCF, so that
the files open in a genome browser and work with bcftools and bedtools.

How the canSNP is carried:

- ``REF`` is the reference base and ``ALT`` the other allele, as the
  specification requires.
- ``INFO/AA`` is the ancestral allele, a key the specification reserves for
  exactly this. The derived allele is the other one, so polarity is explicit
  whether or not the reference carries the derived base.
- ``INFO/NODE`` is the node of the scheme whose genomes all carry the derived
  allele. A SNP that marks no node has neither key.
- The scores are INFO fields, and a canSNP that does not pass the cutoff, is
  not among the best of its node, or lies in a masked region is kept with a
  FILTER code saying so rather than dropped.
- ``##contig`` lines give each contig's ID and length, so that a file is
  checked against the reference it is used with.

Only what this tool writes is read back: biallelic SNPs, one sample column or
none. Reserved characters in INFO values are percent-encoded, as VCF 4.3 asks.
"""

import logging
import os
import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from datetime import date

from cansnpscorer import __version__
from cansnpscorer.model import Reference

VERSION = "VCFv4.3"
COLUMNS = ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
BASES = frozenset("ACGT")

#: The characters VCF 4.3 asks to be percent-encoded in INFO values.
RESERVED = {":": "%3A", ";": "%3B", "=": "%3D", "%": "%25", ",": "%2C",
            "\r": "%0D", "\n": "%0A", "\t": "%09", " ": "%20"}

#: A contig ID as VCF 4.3 allows it (section 1.4.7).
CONTIG_ID = re.compile(r"^[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*$")


def encode(value: str) -> str:
    return "".join(RESERVED.get(character, character) for character in value)


def decode(value: str) -> str:
    if "%" not in value:
        return value
    out = []
    index = 0
    while index < len(value):
        digits = value[index + 1:index + 3]
        if value[index] == "%" and len(digits) == 2 and \
                all(c in "0123456789ABCDEFabcdef" for c in digits):
            out.append(chr(int(value[index + 1:index + 3], 16)))
            index += 3
        else:
            out.append(value[index])
            index += 1
    return "".join(out)


class VcfError(ValueError):
    """The file is not VCF, or holds something this tool did not write."""


@dataclass(frozen=True)
class Info:
    """An INFO field as the header declares it."""

    id: str
    number: str
    type: str
    description: str

    def header(self) -> str:
        description = self.description.replace("\\", "\\\\").replace('"', '\\"')
        return (f'##INFO=<ID={self.id},Number={self.number},Type={self.type},'
                f'Description="{description}">')


@dataclass(frozen=True)
class Filter:
    id: str
    description: str

    def header(self) -> str:
        description = self.description.replace("\\", "\\\\").replace('"', '\\"')
        return f'##FILTER=<ID={self.id},Description="{description}">'


ANCESTRAL = Info("AA", "1", "String", "Ancestral allele")
NODE = Info("NODE", "1", "String",
            "Node of the scheme whose genomes all carry the derived allele, the allele that is "
            "not AA, and that no other genome carries")


#: Distances are reported up to this many bases; what lies further is not known.
NEAREST_CAP = 1000
NEAREST = Info("NEAREST", "1", "Integer",
               f"Distance in bases to the nearest other variation among the genomes placed in "
               f"the scheme: a column where two of them carry different bases, a third allele "
               f"or a partly missing column included, or either end of an insertion or deletion "
               f"among them. A column where the reference has a gap counts at the reference "
               f"bases on either side of it. Up to {NEAREST_CAP}, which means at least "
               f"{NEAREST_CAP}")

MISSING = Info("MISSING", "1", "Integer",
               "Number of genomes of the alignment without a base (gap or N) at the SNP; "
               "0 unless extract was run with --min-called below 1")

COPIES = Info("COPIES", "1", "Integer",
              "Highest number of times a k-mer holding the SNP occurs in the reference, both "
              "strands counted; above 1 the SNP lies in repeated sequence and is ranked after "
              "those in unique sequence. 0 when no k-mer of A, C, G and T holds it")


LOCUS = Info("LOCUS", "1", "String",
             "First and last base on the contig of the gene the SNP lies in; canSNPs of a "
             "node are chosen from different genes where they can be")

OBS_DISCORDANT = Info("OBS_DISCORDANT", "1", "Integer",
                      "Number of genomes of a typing report in which the SNP disagreed with "
                      "the other canSNPs of its node or with the genome's path. At a node "
                      "whose canSNPs are derived in groups of genomes of different width, none "
                      "of them holding more than half of the canSNPs: the number of genomes by "
                      "which the SNP differs from the narrowest group")
OBS_GENOMES = Info("OBS_GENOMES", "1", "Integer",
                   "Number of genomes of the typing report the SNP was observed in; 0 when "
                   "the scheme the report typed with did not hold it")
OBS_MISSING = Info("OBS_MISSING", "1", "Integer",
                   "Number of genomes of a typing report without a base at the SNP")


def count_of(record: "Record", key: str) -> int:
    value = record.info.get(key)
    return int(value) if value not in (None, "") else 0


def copies_of(record: "Record") -> int:
    value = record.info.get("COPIES")
    return int(value) if value not in (None, "") else 1


@dataclass
class Record:
    """One SNP."""

    chrom: str
    pos: int
    """1-based position on the contig."""
    ref: str
    alt: str
    info: dict[str, str | None] = field(default_factory=dict)
    filters: list[str] = field(default_factory=list)
    """The FILTER codes it failed; empty when it passed."""

    @property
    def node(self) -> str | None:
        return self.info.get("NODE")

    @property
    def ancestral(self) -> str | None:
        return self.info.get("AA")

    @property
    def derived(self) -> str | None:
        if self.ancestral is None:
            return None
        return self.alt if self.ancestral == self.ref else self.ref

    @property
    def passed(self) -> bool:
        return not self.filters

    def line(self) -> str:
        fields = []
        for key, value in self.info.items():
            fields.append(key if value is None else f"{key}={encode(str(value))}")
        return "\t".join([
            self.chrom, str(self.pos), ".", self.ref, self.alt, ".",
            ";".join(self.filters) if self.filters else "PASS",
            ";".join(fields) if fields else ".",
        ])


def check_contig_names(reference: Reference) -> None:
    """Refuse contig names VCF 4.3 does not allow as IDs, such as one holding a
    comma, a bracket or whitespace. Called before any work is done, so that a
    run does not fail only when its output is written.

    Raises:
        VcfError: If any contig name is not allowed.
    """
    illegal = [contig.name for contig in reference if not CONTIG_ID.match(contig.name)]
    if illegal:
        raise VcfError(
            f"the reference names {len(illegal)} contigs VCF 4.3 does not allow as IDs: "
            f"{', '.join(repr(name) for name in illegal[:3])}. Rename them in the reference "
            f"FASTA, which the alignment and annotations must then follow")


def header_lines(reference: Reference, reference_path: str, command: str,
                 infos: Sequence[Info], filters: Sequence[Filter] = (),
                 extra: Sequence[str] = ()) -> list[str]:
    lines = [
        f"##fileformat={VERSION}",
        f"##fileDate={date.today().strftime('%Y%m%d')}",
        f"##source=CanSNPCreator {__version__} {command}",
        f"##reference=file://{os.path.abspath(reference_path)}",
    ]
    check_contig_names(reference)
    lines += [f"##contig=<ID={contig.name},length={contig.length}>" for contig in reference]
    lines += [info.header() for info in infos]
    lines += [definition.header() for definition in filters]
    lines += list(extra)
    lines.append("\t".join(COLUMNS))
    return lines


def write(path: str, header: Sequence[str], records: Iterable[Record]) -> int:
    """Write a VCF file. Returns the number of records written."""
    count = 0
    with open(path, "w", newline="\n") as handle:
        handle.writelines(line + "\n" for line in header)
        for record in records:
            handle.write(record.line() + "\n")
            count += 1
    return count


@dataclass
class Vcf:
    meta: list[str]
    contigs: dict[str, int]
    infos: dict[str, Info]
    records: list[Record]

    def check_against(self, reference: Reference, path: str) -> None:
        """The file's contigs must be the reference's, by ID and length."""
        if not self.contigs:
            raise VcfError(f"{path}: declares no ##contig lines, so it cannot be checked "
                           f"against the reference")
        expected = {contig.name: contig.length for contig in reference}
        if self.contigs != expected:
            differ = [
                f"{name}: {self.contigs.get(name, 'absent')} in the VCF, "
                f"{expected.get(name, 'absent')} in the reference"
                for name in sorted(set(self.contigs) | set(expected))
                if self.contigs.get(name) != expected.get(name)
            ]
            raise VcfError(
                f"{path}: was written against another reference; contig lengths differ for "
                + "; ".join(differ[:3]))
        wrong = []
        for record in self.records:
            if record.chrom not in expected:
                wrong.append(f"{record.chrom} is not a contig of the reference")
            elif not 1 <= record.pos <= expected[record.chrom]:
                wrong.append(f"{record.chrom}:{record.pos} is past the end of the contig")
            else:
                base = reference.contig(record.chrom).base(record.pos - 1).upper()
                if base != record.ref:
                    wrong.append(f"{record.chrom}:{record.pos} has REF {record.ref} where the "
                                 f"reference has {base}")
        if wrong:
            raise VcfError(
                f"{path}: {len(wrong)} records do not match the reference: "
                + "; ".join(wrong[:3]))


def _attributes(text: str) -> dict[str, str]:
    """The key=value pairs of a structured header line, quotes removed."""
    attributes: dict[str, str] = {}
    index = 0
    while index < len(text):
        equals = text.index("=", index)
        key = text[index:equals]
        index = equals + 1
        if index < len(text) and text[index] == '"':
            index += 1
            value = []
            while text[index] != '"':
                if text[index] == "\\":
                    index += 1
                value.append(text[index])
                index += 1
            index += 1
            attributes[key] = "".join(value)
        else:
            end = text.find(",", index)
            end = len(text) if end < 0 else end
            attributes[key] = text[index:end]
            index = end
        if index < len(text) and text[index] == ",":
            index += 1
    return attributes


def read(path: str) -> Vcf:
    """Read a VCF file of biallelic SNPs.

    Raises:
        VcfError: If the file is not VCF, declares no INFO it uses, or holds a
            record that is not a biallelic SNP.
        OSError: If the file cannot be read.
    """
    meta: list[str] = []
    contigs: dict[str, int] = {}
    infos: dict[str, Info] = {}
    records: list[Record] = []
    seen_header = False

    from cansnpscorer.helpers.file_utils import open_text

    with open_text(path, newline="") as handle:
        for number, line in enumerate(handle, 1):
            line = line.rstrip("\r\n")
            if number == 1 and not line.startswith("##fileformat=VCF"):
                raise VcfError(f"{path}: is not VCF; its first line is not ##fileformat")
            if line.startswith("##"):
                meta.append(line)
                try:
                    if line.startswith("##contig=<") and line.endswith(">"):
                        attributes = _attributes(line[len("##contig=<"):-1])
                        if "ID" in attributes and "length" in attributes:
                            contigs[attributes["ID"]] = int(attributes["length"])
                    elif line.startswith("##INFO=<") and line.endswith(">"):
                        attributes = _attributes(line[len("##INFO=<"):-1])
                        infos[attributes["ID"]] = Info(
                            attributes["ID"], attributes.get("Number", "."),
                            attributes.get("Type", "String"), attributes.get("Description", ""))
                except (KeyError, IndexError, ValueError):
                    raise VcfError(f"{path}: line {number} is a malformed header line") from None
                continue
            if line.startswith("#CHROM"):
                if tuple(line.split("\t")[:8]) != COLUMNS:
                    raise VcfError(f"{path}: line {number} is not the VCF column header")
                seen_header = True
                continue
            if not line.strip():
                continue
            if not seen_header:
                raise VcfError(f"{path}: line {number} comes before the #CHROM header")

            fields = line.split("\t")
            if len(fields) < 8:
                raise VcfError(f"{path}: line {number} holds {len(fields)} of the 8 columns")
            chrom, pos, _id, ref, alt, _qual, filters, info = fields[:8]
            ref, alt = ref.upper(), alt.upper()
            if ref not in BASES or alt not in BASES or ref == alt:
                raise VcfError(
                    f"{path}: line {number} is not a biallelic SNP ({ref} to {alt}); "
                    f"only SNPs as extract writes them are read")
            values: dict[str, str | None] = {}
            if info != ".":
                for item in info.split(";"):
                    key, has_value, value = item.partition("=")
                    if key not in infos:
                        raise VcfError(
                            f"{path}: line {number} uses INFO/{key}, which the header does not "
                            f"declare")
                    values[key] = decode(value) if has_value else None
            try:
                position = int(pos)
            except ValueError:
                raise VcfError(f"{path}: line {number} has the position '{pos}'") from None
            failed = [] if filters in ("PASS", ".") else filters.split(";")
            records.append(Record(chrom, position, ref, alt, values, failed))

    if not seen_header:
        raise VcfError(f"{path}: has no #CHROM header line")
    logging.info(f"Read {len(records)} SNPs from {path}")
    return Vcf(meta, contigs, infos, records)


def snp_record(reference: Reference, contig: str, position: int, alleles: tuple[str, str],
               node: str | None = None, ancestral: str | None = None) -> Record:
    """A record for a SNP at a 1-based position, REF taken from the reference.

    Raises:
        VcfError: If the reference base is neither allele.
    """
    ref = reference.contig(contig).base(position - 1).upper()
    if ref not in alleles:
        raise VcfError(
            f"{contig}:{position} is {ref} in the reference, which is neither allele "
            f"{alleles[0]} nor {alleles[1]}")
    alt = alleles[1] if alleles[0] == ref else alleles[0]
    info: dict[str, str | None] = {}
    if node is not None:
        info["NODE"] = node
        info["AA"] = ancestral
    return Record(contig, position, ref, alt, info)


# -- the files of extract and evaluate ---------------------------------------

SCORE = Info("SCORE", "1", "Integer",
             "Weighted sum of the five score terms, weighted as ##CanSNPCreatorWeights gives; "
             "higher is a better context for calling the SNP")


def term_infos() -> list[Info]:
    from cansnpscorer.helpers.scoring_methods import INFO_KEYS

    descriptions = {
        "neighbour_score": "distance to the closest neighbouring SNP",
        "gc_score": "closeness of the surrounding GC content to the goal",
        "tandem_score": "distance to the closest tandem repeat",
        "gene_score": "whether the SNP lies in a gene or CDS",
        "mutation_score": "100 if synonymous in its CDS, 50 if it changes the amino acid, 0 if "
                          "it gains or loses a stop codon or lies outside every CDS",
    }
    return [Info(key, "1", "Integer", f"Score term, 0 to 100 times its weight: "
                                      f"{descriptions[term]}")
            for term, key in INFO_KEYS.items()]


def score_filters(cutoff: int | None, top: int | None, spacing: int = 0) -> list[Filter]:
    filters = [Filter("LowScore", f"SCORE is not above the cutoff of {cutoff}"
                                  if cutoff is not None else "Not applied: no cutoff was set")]
    if top is not None:
        apart = (f", and none within {spacing} bases of one already chosen"
                 if spacing > 0 else "")
        filters.append(Filter(
            "NotTop", f"Neither LowScore nor Masked, but not among the {top} canSNPs chosen for "
                      f"its node. These are taken by reliability first (copy number in the "
                      f"reference, disagreement in a typing report, whether a typing report "
                      f"observed it, missing bases), then by SCORE, from different genes where "
                      f"they can be{apart}"))
    filters.append(Filter("Masked", "Lies in a region given as masked, such as a recombinant "
                                    "region"))
    return filters


@dataclass
class Snp:
    """A SNP as the extract subcommand writes it."""

    contig: str
    position: int
    node: str | None
    derived: str
    ancestral: str
    nearest: int | None = None
    """Distance to the nearest variation among the scheme's genomes, when known."""
    missing: int = 0
    """Genomes of the alignment without a base at the SNP."""


def nearest_of(record: Record) -> int | None:
    value = record.info.get("NEAREST")
    return int(value) if value not in (None, "") else None


def missing_of(record: Record) -> int:
    value = record.info.get("MISSING")
    return int(value) if value not in (None, "") else 0


def read_snps(path: str, reference: Reference) -> list[Snp]:
    """Read the SNPs extract wrote, checked against the reference.

    For a SNP that marks no node the polarity is unknown, and the REF and ALT
    alleles are given as derived and ancestral only to carry both bases.
    """
    parsed = read(path)
    parsed.check_against(reference, path)
    snps = []
    for record in parsed.records:
        if record.node is not None:
            if record.ancestral not in (record.ref, record.alt):
                raise VcfError(
                    f"{path}: {record.chrom}:{record.pos} gives the ancestral allele "
                    f"{record.ancestral}, which is neither REF nor ALT")
            snps.append(Snp(record.chrom, record.pos, record.node, record.derived,
                            record.ancestral, nearest_of(record), missing_of(record)))
        else:
            snps.append(Snp(record.chrom, record.pos, None, record.alt, record.ref,
                            nearest_of(record), missing_of(record)))
    if not snps:
        raise VcfError(f"{path}: holds no SNPs")
    return snps


def write_scores(path: str, scores, reference: Reference, reference_path: str,
                 cutoff: int | None, top: int | None, weights: dict[str, float],
                 copies: bool = False, observed: bool = False, spacing: int = 0) -> int:
    """Write scored canSNPs, each with its terms and the filters it failed;
    with ``copies``, with the copy number of their context as well, and with
    ``observed``, with how they behaved in a typing report. ``spacing`` is
    the fewest bases between two canSNPs chosen for a node, which the
    description of NotTop names."""
    from cansnpscorer.helpers.scoring_methods import INFO_KEYS, TERMS

    records = []
    for score in sorted(scores, key=lambda s: (reference.contig(s.contig).start, s.position)):
        record = snp_record(reference, score.contig, score.position,
                            (score.derived, score.ancestral), score.node, score.ancestral)
        record.info["SCORE"] = str(score.total)
        for term in TERMS:
            record.info[INFO_KEYS[term]] = str(score.terms[term])
        record.info["MISSING"] = str(getattr(score, "missing", 0))
        if copies:
            record.info["COPIES"] = str(getattr(score, "copies", 1))
        if getattr(score, "locus", ""):
            record.info["LOCUS"] = score.locus
        if observed:
            record.info["OBS_GENOMES"] = str(getattr(score, "observed_genomes", 0))
            record.info["OBS_DISCORDANT"] = str(getattr(score, "discordant", 0))
            record.info["OBS_MISSING"] = str(getattr(score, "observed_missing", 0))
        record.filters = list(score.filters)
        records.append(record)

    weighting = ",".join(f"{term}={weights.get(term, 1.0):g}" for term in TERMS)
    header = header_lines(reference, reference_path, "evaluate",
                          [NODE, ANCESTRAL, SCORE, *term_infos(), MISSING,
                           *([COPIES] if copies else []), LOCUS,
                           *([OBS_GENOMES, OBS_DISCORDANT, OBS_MISSING] if observed else [])],
                          score_filters(cutoff, top, spacing),
                          [f"##CanSNPCreatorWeights={weighting}"])
    return write(path, header, records)


def read_scores(path: str, reference: Reference | None = None):
    """Read the scored canSNPs evaluate wrote."""
    from cansnpscorer.helpers.scoring_methods import INFO_KEYS, ScoredSnp

    parsed = read(path)
    if reference is not None:
        parsed.check_against(reference, path)
    missing = [key for key in ("NODE", "AA", "SCORE", *INFO_KEYS.values())
               if key not in parsed.infos]
    if missing:
        raise VcfError(
            f"{path}: declares no INFO/{', INFO/'.join(missing)}; it should be the output of "
            f"the evaluate subcommand")
    scores = []
    for record in parsed.records:
        if record.node is None:
            continue
        if record.ancestral not in (record.ref, record.alt):
            raise VcfError(
                f"{path}: {record.chrom}:{record.pos} gives the ancestral allele "
                f"{record.ancestral}, which is neither REF nor ALT")
        try:
            scores.append(ScoredSnp(
                contig=record.chrom, position=record.pos, node=record.node,
                total=int(record.info["SCORE"]),
                terms={term: int(record.info[key]) for term, key in INFO_KEYS.items()},
                derived=record.derived, ancestral=record.ancestral,
                filters=list(record.filters), missing=missing_of(record),
                copies=copies_of(record), locus=record.info.get("LOCUS") or "",
                discordant=count_of(record, "OBS_DISCORDANT"),
                observed_genomes=count_of(record, "OBS_GENOMES"),
                observed_missing=count_of(record, "OBS_MISSING")))
        except (KeyError, TypeError, ValueError):
            raise VcfError(
                f"{path}: {record.chrom}:{record.pos} lacks a whole number score") from None
    if not scores:
        raise VcfError(f"{path}: holds no scored canSNPs")
    logging.info(f"Read {len(scores)} scored canSNPs, {sum(s.passed for s in scores)} passing, "
                 f"from {path}")
    return scores
