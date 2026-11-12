"""Checks on the inputs, reported together rather than one failure at a time.

The readers refuse what they cannot read at all. The checks here are the ones
that need more than one file to see: whether the alignment, the group file, the
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
from cansnpscorer.model import Reference, ungap

Report = Tuple[List[str], List[str]]


def check_alignment_against_reference(alignment: Mapping[str, str],
                                     reference: Reference) -> Report:
    """The alignment's first record must be this reference genome."""
    problems: List[str] = []
    warnings: List[str] = []

    if not alignment:
        return ["the alignment holds no genomes"], warnings

    name, aligned = next(iter(alignment.items()))
    ungapped = ungap(aligned)

    if len(ungapped) != reference.length:
        problems.append(
            f"the alignment's first record '{name}' holds {len(ungapped)} bases once its gaps "
            f"are removed, but the reference genome holds {reference.length}. The first record "
            f"of the alignment must be the reference genome, and the reference file must hold "
            f"the contigs the alignment concatenates, in that order")
    elif ungapped.upper() != reference.sequence.upper():
        problems.append(
            f"the alignment's first record '{name}' is the same length as the reference genome "
            f"but not the same sequence")

    gaps = len(aligned) - len(ungapped)
    if gaps:
        warnings.append(
            f"the alignment's reference row carries {gaps} gaps, so {gaps} columns have no "
            f"position in the reference genome and any SNP in them is left out")

    return problems, warnings


def check_groups(alignment: Mapping[str, str], clusters: Mapping[str, Sequence[str]],
                 hierarchy: Hierarchy, reference_name: str = "") -> Report:
    """The group file, the hierarchy and the alignment must describe one set of genomes."""
    problems: List[str] = []
    warnings: List[str] = []

    grouped: Dict[str, List[str]] = {}
    for cluster, isolates in clusters.items():
        for isolate in isolates:
            grouped.setdefault(isolate, []).append(cluster)

    missing = sorted(set(grouped) - set(alignment))
    if missing:
        problems.append(
            f"the group file names {len(missing)} isolates that the alignment does not hold: "
            f"{', '.join(missing[:5])}")

    unknown = sorted(set(clusters) - set(hierarchy.parents))
    if unknown:
        problems.append(
            f"the group file names {len(unknown)} clusters that the hierarchy does not hold: "
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
            f"{len(ungrouped)} genomes of the alignment are in no cluster, and they still count "
            f"against every canSNP: {', '.join(ungrouped[:5])}")

    known = {name: isolates for name, isolates in clusters.items() if name in hierarchy}
    without_isolates = sorted(
        node for node, isolates in hierarchy.isolates(known).items() if not isolates)
    if without_isolates:
        warnings.append(
            f"{len(without_isolates)} nodes of the hierarchy have no isolates, so no canSNP can "
            f"mark them: {', '.join(without_isolates[:5])}")

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

    empty = sorted(name for name in reference.names if not annotations.get(name))
    if empty and annotations:
        warnings.append(
            f"these contigs have no {description} annotation, so every SNP on them scores zero "
            f"for it: {', '.join(empty[:5])}")

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
        f"The {description} are consistent"
        + (f", with {len(warnings)} warnings" if warnings else ""))
    return True
