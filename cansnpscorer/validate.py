"""Checks on the inputs, reported together rather than one failure at a time.

The readers refuse what they cannot read at all. The checks here are the ones
that need more than one file to see: whether the alignment, the genome map, the
hierarchy and the reference describe the same genomes and the same genome. They
collect problems and warnings in lists, in the same shape as the validation
script of the CanSNPer2-data repository, so that a run reports everything that
is wrong with its inputs before it starts work rather than stopping at the
first thing it meets.

A problem makes the result wrong or impossible; a warning leaves the run usable
but says something about the data worth knowing, most often a part of the
scheme that cannot be typed.
"""

import logging
from typing import Dict, List, Mapping, Sequence, Tuple

from cansnpscorer.hierarchy import Hierarchy
from cansnpscorer.model import AlignmentMap, Reference, ReferenceError, ungap

Report = Tuple[List[str], List[str]]


def check_alignment_against_reference(alignment: Mapping[str, str],
                                     reference: Reference) -> Report:
    """The alignment's first record must be this reference genome."""
    problems: List[str] = []
    warnings: List[str] = []

    if not alignment:
        return ["the alignment holds no genomes"], warnings

    name, aligned = next(iter(alignment.items()))

    try:
        mapping = AlignmentMap(aligned, reference)
    except ReferenceError as error:
        problems.append(f"the alignment's first record '{name}' is not the reference: {error}")
        return problems, warnings

    if mapping.projected:
        if mapping.unaligned:
            warnings.append(
                f"the alignment is projected onto reference coordinates, and "
                f"{mapping.unaligned} positions of the reference are covered by no aligned "
                f"block, so no SNP can be found there")
    else:
        gaps = len(aligned) - len(ungap(aligned))
        if gaps:
            warnings.append(
                f"the alignment's reference row carries {gaps} gaps, so {gaps} columns have no "
                f"position in the reference genome and any SNP in them is left out")

    return problems, warnings


def check_groups(alignment: Mapping[str, str], clusters: Mapping[str, Sequence[str]],
                 hierarchy: Hierarchy, reference_name: str = "") -> Report:
    """The genome map, the hierarchy and the alignment must describe one set of genomes."""
    problems: List[str] = []
    warnings: List[str] = []

    grouped: Dict[str, List[str]] = {}
    for cluster, isolates in clusters.items():
        for isolate in isolates:
            grouped.setdefault(isolate, []).append(cluster)

    missing = sorted(set(grouped) - set(alignment))
    if missing:
        problems.append(
            f"the genome map names {len(missing)} isolates that the alignment does not hold: "
            f"{', '.join(missing[:5])}")

    unknown = sorted(set(clusters) - set(hierarchy.parents))
    if unknown:
        problems.append(
            f"the genome map names {len(unknown)} clusters that the hierarchy does not hold: "
            f"{', '.join(unknown[:5])}")

    repeated = sorted(name for name, found in grouped.items() if len(found) > 1)
    if repeated:
        examples = "; ".join(
            f"{name} in {' and '.join(grouped[name])}" for name in repeated[:3])
        problems.append(
            f"{len(repeated)} isolates are in more than one cluster, so no count of them can be "
            f"right: {examples}")

    # The reference genome is normally not a member of the scheme, so it is
    # expected to be ungrouped; anything else that is ungrouped silently
    # weakens every canSNP, because a genome outside the scheme still carries
    # bases that decide whether a SNP is canonical.
    ungrouped = sorted(set(alignment) - set(grouped) - {reference_name})
    if ungrouped:
        warnings.append(
            f"{len(ungrouped)} genomes of the alignment are placed at no node, so they count as "
            f"outside every node, as an outgroup should; any other genome among them weakens "
            f"the canSNPs it differs at: {', '.join(ungrouped[:5])}")

    known = {name: isolates for name, isolates in clusters.items() if name in hierarchy}
    without_isolates = sorted(
        node for node, isolates in hierarchy.isolates(known).items() if not isolates)
    if without_isolates:
        warnings.append(
            f"{len(without_isolates)} nodes of the hierarchy have no isolates, so no canSNP can "
            f"mark them: {', '.join(without_isolates[:5])}")

    return problems, warnings


def unresolved_genomes(unmatched: Sequence[str], ambiguous: Sequence[str]) -> Report:
    """Genome names that match no alignment record, or several.

    A genome the alignment does not hold is left out with a warning: a genome
    map from RepGenR2 names every genome, while its alignment often holds only
    the representatives. A name that matches two records is a problem, since
    either choice could be wrong.
    """
    problems: List[str] = []
    warnings: List[str] = []
    if unmatched:
        warnings.append(
            f"the genome map names {len(unmatched)} genomes that the alignment does not hold, "
            f"and they were left out: {', '.join(unmatched[:5])}")
    if ambiguous:
        problems.append(
            f"{len(ambiguous)} genomes of the genome map match more than one record of the "
            f"alignment: {', '.join(ambiguous[:5])}")
    return problems, warnings


def check_annotations(annotations: Mapping[str, Sequence[Sequence[int]]], reference: Reference,
                      description: str) -> Report:
    """Annotations must lie within the contigs they are attached to."""
    problems: List[str] = []
    warnings: List[str] = []

    for name, intervals in annotations.items():
        if name not in reference.names:
            continue
        contig = reference.contig(name)
        outside = [
            interval for interval in intervals
            if interval[0] < 0 or interval[1] >= contig.length
        ]
        if outside:
            problems.append(
                f"{len(outside)} {description} intervals on '{name}' lie outside the contig, "
                f"which is {contig.length} bases long; the first is "
                f"{outside[0][0] + 1}..{outside[0][1] + 1}")

    # A contig without genes is worth knowing, since its SNPs score nothing for
    # the gene and mutation terms. A contig without repeats is ordinary, and
    # its SNPs score the most for the tandem term.
    empty = sorted(name for name in reference.names if not annotations.get(name))
    if empty and annotations and description == "gene":
        warnings.append(
            f"these contigs have no gene annotation, so every SNP on them scores zero for the "
            f"gene and mutation terms: {', '.join(empty[:5])}")

    return problems, warnings


def report(problems: Sequence[str], warnings: Sequence[str], description: str = "inputs") -> bool:
    """Log a report. Returns True when the run may go ahead."""
    for warning in warnings:
        logging.warning(warning)
    for problem in problems:
        logging.error(problem)

    if problems:
        logging.error(
            f"{len(problems)} problems with the {description}; nothing was written")
        return False

    logging.info(
        f"No problems found in the {description}"
        + (f", with {len(warnings)} warnings" if warnings else ""))
    return True
