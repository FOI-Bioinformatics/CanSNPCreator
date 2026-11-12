"""Cleaning of sequence and strain names.

Names reach the tools from FASTA headers, group files and annotation files, and
the same genome is often written differently in each. They are reduced to their
alphanumeric characters so that the three can be matched, which is also why a
cluster name must be compared before it is cleaned: a cluster called A.1 must
stay A.1 to match the hierarchy.
"""

import logging
import re

NOT_ALPHANUMERIC = re.compile("[^A-Za-z0-9]+")
FASTA_SUFFIX = re.compile(r"\.(fasta|fna|fa|fas|ffn)(\.gz)?$", re.IGNORECASE)


#: An NCBI assembly accession once cleaned, as in GCA0000071251.
ACCESSION = re.compile(r"GC[AF]\d{9}\d+")


def resolve_names(names, records):
    """Match genome names from one file to the record names of another.

    A name matches a record when the two are equal once cleaned. Failing
    that, a name that is an assembly accession matches the one record whose
    name contains it, since RepGenR2 names records ``Family_Genus_species_Acc``
    while its genomes_map.tsv gives the accession alone.

    Args:
        names: The names to match, as the file gives them.
        records: The record names to match them against, already cleaned.

    Returns:
        A dictionary from each name that matched to its record, the names that
        matched nothing, and the names that matched more than one record.
    """
    by_clean = {clean_name(record): record for record in records}
    matched = {}
    unmatched = []
    ambiguous = []
    for name in names:
        cleaned = clean_name(strip_fasta_suffix(name))
        if cleaned in by_clean:
            matched[name] = by_clean[cleaned]
            continue
        if ACCESSION.fullmatch(cleaned):
            # Not followed by a digit, so that version 1 does not match version 10.
            inside = re.compile(re.escape(cleaned) + r"(?!\d)")
            found = [record for key, record in by_clean.items() if inside.search(key)]
            if len(found) == 1:
                matched[name] = found[0]
                continue
            if found:
                ambiguous.append(name)
                continue
        unmatched.append(name)
    return matched, unmatched, ambiguous


def resolve_clusters(clusters, records):
    """Rewrite the genomes of each node as the alignment records they match.

    Returns the clusters with record names, the genome names that matched no
    record, and those that matched more than one.
    """
    names = [name for genomes in clusters.values() for name in genomes]
    matched, unmatched, ambiguous = resolve_names(names, records)
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
    Reduce a name to its alphanumeric characters.

    Args:
        name: The name as the file gives it.

    Returns:
        The name with every other character removed.

    Example:
        >>> clean_name("GCA_000007125.1 ASM712v1")
        'GCA0000071251ASM712v1'
    """
    cleaned = NOT_ALPHANUMERIC.sub("", name)
    if cleaned != name:
        logging.debug(f"Cleaned the name {name!r} to {cleaned!r}")
    return cleaned
