"""The collapse and place subcommands. The work is in tree_tools."""

import logging

from cansnpscorer.findings import PROBLEM
from cansnpscorer.helpers.file_utils import read_alignment, reference_parser
from cansnpscorer.hierarchy import HierarchyError
from cansnpscorer.model import ReferenceError
from cansnpscorer.scheme import SchemeError, read_scheme
from cansnpscorer.tree_tools import (
    collapse,
    place,
    read_genome_tree,
    read_members,
    rename,
    write_collapsed,
    write_placement,
)


def collapse_main(tree_file=None, nodes_file=None, members=None, min_support=None,
                  min_genomes=2, min_length=None, prefix_names=None, outdir=".",
                  prefix="collapsed", outgroups=()) -> int:
    """Collapse a genome tree into a scheme's hierarchy. Returns an exit status."""
    try:
        tree = read_genome_tree(tree_file, nodes_file)
        folded = read_members(members, tree) if members else None
        collapsed = collapse(tree, min_support, min_genomes, min_length, folded, outgroups)
        if prefix_names:
            collapsed = rename(collapsed, prefix_names)
        written = write_collapsed(collapsed, tree, outdir, prefix)
    except (OSError, ValueError, HierarchyError) as error:
        logging.error(error)
        return 1

    if min_support is not None and not tree.support:
        logging.warning(
            "--min-support was given but the tree carries no support values, so no node was "
            "removed for its support")
    def kind(part: str) -> str:
        if part.startswith("support"):
            return "support"
        if part.startswith("branch length"):
            return "branch length"
        if part.endswith(" genomes"):
            return "too few genomes"
        return "the same genomes as the parent"

    reasons = {}
    for reason in collapsed.removed.values():
        for part in reason.split(", "):
            reasons[kind(part)] = reasons.get(kind(part), 0) + 1
    logging.info(
        f"Removed {len(collapsed.removed)} internal nodes "
        f"({', '.join(f'{count} for {key}' for key, count in sorted(reasons.items()))}); "
        f"wrote {', '.join(written)}")
    return 0


def place_main(scheme_dir=None, organism=None, tree_file=None, alignment=None,
               reference=None, genome=None, members=None, min_match=0.5, output=None,
               report=None) -> int:
    """Place the genomes of a tree on an existing scheme. Returns an exit status."""
    try:
        scheme = read_scheme(scheme_dir, organism)
        tree = read_genome_tree(tree_file)
        folded = read_members(members, tree) if members else None
        aligned = read_alignment(alignment)
        sequence = reference_parser(reference)
        placement = place(scheme, tree, aligned, sequence, genome, min_match, folded)
        write_placement(placement, output, report)
    except (OSError, ValueError, HierarchyError, SchemeError, ReferenceError) as error:
        logging.error(error)
        return 1

    for finding in placement.findings:
        (logging.error if finding.severity == PROBLEM else logging.warning)(str(finding))
    logging.info(f"Wrote the genome map to {output}" + (f" and the report to {report}"
                                                         if report else ""))
    return 2 if any(finding.severity == PROBLEM for finding in placement.findings) else 0
