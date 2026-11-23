"""Matching genome names written differently in different files.

Sequences are named by their FASTA ID, exactly, as GFF3, VCF and BED name
them; nothing here renames a sequence. A genome, though, is often written one
way in an alignment, another in a genome map and a third in a tree: a file
name, a record ID, or an assembly accession alone. ``resolve_names`` matches
such names to alignment records once their separators are written alike, or
by an accession the record name holds, and the record ID is what is used from
then on. Node names are never cleaned: a node called A.1 stays A.1.
"""

import logging
import re

NOT_ALPHANUMERIC = re.compile("[^A-Za-z0-9]+")
FASTA_SUFFIX = re.compile(r"\.(fasta|fna|fa|fas|ffn)(\.gz)?$", re.IGNORECASE)


#: An NCBI assembly accession, with or without its version.
ACCESSION = re.compile(r"(GC[AF]_\d{9})(?:\.(\d+))?")
#: An accession inside a longer name, delimited on both sides: not part of a
#: longer word before it, and not followed by a further digit. Its version is
#: every digit after the dot, so that version 1 is not version 10.
HELD_ACCESSION = re.compile(r"(?<![A-Za-z0-9])(GC[AF]_\d{9})(?!\d)(?:\.(\d+))?")


def resolve_names(names, records):
    """Match genome names from one file to the record names of another.

    In order: a name equal to a record matches it, then one equal to it once
    a FASTA suffix is removed from either. Failing that, a name matches the
    record it equals once every run of characters other than letters and
    digits is written as one separator, provided only one record does; the
    separators are kept, so that ``S1_2`` is not ``S12``. Failing that, a name
    that is an assembly accession matches the one record whose name contains
    it, since RepGenR2 names records ``Family_Genus_species_Acc`` while its
    genomes_map.tsv gives the accession alone. The version must be the same,
    so version 1 is not version 10, and what follows it in the record name
    does not matter, as in an NCBI file name ``GCA_900475405.1_49699_E01``. An
    accession without a version matches the record that holds it in any
    version, when only one does. A name that would match two records is
    reported, not guessed.

    Args:
        names: The names to match, as the file gives them.
        records: The record names to match them against, as their IDs.

    Returns:
        A dictionary from each name that matched to its record, the names that
        matched nothing, and the names that matched more than one record.
    """
    by_name = {}
    by_clean = {}
    for record in records:
        for key in {record, strip_fasta_suffix(record)}:
            by_name.setdefault(key, set()).add(record)
        by_clean.setdefault(clean_name(strip_fasta_suffix(record)), set()).add(record)

    # The accessions each record name holds, delimited as described above, so
    # that an accession is looked up rather than searched for in every record.
    by_stem = {}
    by_version = {}
    for record in records:
        for found in HELD_ACCESSION.finditer(record):
            stem, version = found.groups()
            by_stem.setdefault(stem, set()).add(record)
            if version is not None:
                by_version.setdefault((stem, version), set()).add(record)

    exact = set(records)
    matched = {}
    unmatched = []
    ambiguous = []
    for name in names:
        stripped = strip_fasta_suffix(name)
        if name in exact:
            # The record's own ID, whatever other record resembles it.
            matched[name] = name
            continue
        if name in by_name and len(by_name[name]) == 1:
            matched[name] = next(iter(by_name[name]))
            continue
        if stripped in by_name and len(by_name[stripped]) == 1:
            matched[name] = next(iter(by_name[stripped]))
            continue
        found = by_clean.get(clean_name(stripped), set())
        accession = ACCESSION.fullmatch(stripped)
        if not found and accession:
            stem, version = accession.groups()
            found = by_version.get((stem, version), set()) if version else by_stem.get(stem, set())
        if len(found) == 1:
            matched[name] = next(iter(found))
        elif found:
            ambiguous.append(name)
        else:
            unmatched.append(name)
    return matched, unmatched, ambiguous


def resolve_clusters(clusters, records):
    """Rewrite the genomes of each node as the alignment records they match.

    Returns the clusters with record names, the genome names that matched no
    record, and those that matched more than one.
    """
    names = [name for genomes in clusters.values() for name in genomes]
    matched, unmatched, ambiguous = resolve_names(names, records)
    # Two names written differently can be one record; placed at two nodes
    # they would put one genome in two places, so both are refused.
    nodes_of: dict[str, set[str]] = {}
    for node, genomes in clusters.items():
        for name in genomes:
            if name in matched:
                nodes_of.setdefault(matched[name], set()).add(node)
    split = {record for record, nodes in nodes_of.items() if len(nodes) > 1}
    for name in list(matched):
        if matched[name] in split:
            ambiguous.append(name)
            del matched[name]
    resolved = {}
    for node, genomes in clusters.items():
        found = []
        for name in genomes:
            record = matched.get(name)
            if record is not None and record not in found:
                found.append(record)
        resolved[node] = found
    return resolved, sorted(set(unmatched)), sorted(set(ambiguous))


def strip_fasta_suffix(name: str) -> str:
    """Remove a FASTA file suffix, so that a file name reads as its genome.

    RepGenR2 names each record and leaf by the genome file name without its
    FASTA suffix, so ``x.fna.gz`` and ``x`` are the same genome.
    """
    return FASTA_SUFFIX.sub("", name)


def clean_name(name: str) -> str:
    """
    Write every run of characters other than letters and digits as one
    underscore, and drop those at either end.

    The separators are kept, not removed: removing them makes ``S1_2`` and
    ``S12`` one name, and one genome would be placed as the other.

    Example:
        >>> clean_name("GCA_000007125.1 ASM712v1")
        'GCA_000007125_1_ASM712v1'
    """
    cleaned = NOT_ALPHANUMERIC.sub("_", name).strip("_")
    if cleaned != name:
        logging.debug(f"Cleaned the name {name!r} to {cleaned!r}")
    return cleaned
