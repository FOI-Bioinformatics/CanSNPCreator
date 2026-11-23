"""The collapse, name and place subcommands. The work is in tree_tools."""

import logging
from dataclasses import dataclass

from cansnpscorer.findings import PROBLEM
from cansnpscorer.helpers.file_utils import (
    fasta_ids,
    read_alignment,
    read_genome_map,
    read_hierarchy,
    read_labels,
    reference_parser,
)
from cansnpscorer.helpers.names import resolve_names
from cansnpscorer.hierarchy import HierarchyError
from cansnpscorer.model import ReferenceError
from cansnpscorer.scheme import Scheme, SchemeError, read_scheme
from cansnpscorer.tree_tools import (
    collapse,
    name_from_labels,
    outgroup_leaves,
    place,
    place_by_clade,
    read_genome_tree,
    read_members,
    rename,
    reroot,
    write_collapsed,
    write_named,
    write_placement,
)


def name_main(hierarchy=None, genome_map=None, labels=None, outdir=".", prefix="named") -> int:
    """Name the nodes of a hierarchy after labelled genomes. Returns an exit status."""
    try:
        tree = read_hierarchy(hierarchy)
        clusters = read_genome_map(genome_map)
        given = read_labels(labels)
    except (OSError, ValueError, HierarchyError) as error:
        logging.error(error)
        return 1

    # A label table usually names genomes by accession and a genome map by
    # record name, so they are matched as every other genome name is.
    placed = [genome for genomes in clusters.values() for genome in genomes]
    matched, unmatched, ambiguous = resolve_names(list(given), placed)
    resolved = {matched[genome]: label for genome, label in given.items() if genome in matched}
    try:
        named = name_from_labels(tree, clusters, resolved)
        written = write_named(named, tree, outdir, prefix)
    except (OSError, ValueError, HierarchyError) as error:
        logging.error(error)
        return 1

    if unmatched or ambiguous:
        logging.warning(
            f"{len(unmatched) + len(ambiguous)} labelled genomes are not in the genome map, or "
            f"match more than one genome of it, and name nothing: "
            f"{', '.join([*unmatched, *ambiguous][:5])}"
        )
    if named.unplaced:
        logging.warning(
            f"{len(named.unplaced)} labels share every node with another label, so no node of "
            f"the scheme tells their genomes apart, and they name nothing: "
            f"{', '.join(named.unplaced)}"
        )
    for label, nodes in sorted(named.split.items()):
        logging.warning(
            f"the genomes labelled {label} do not form one clade, so {len(nodes)} nodes are "
            f"named after it: {', '.join(named.names[node] for node in nodes)}"
        )
    logging.info(
        f"Named {len(named.labels)} of {len(named.names)} nodes after labels; wrote {', '.join(written)}"
    )
    return 0


def collapse_main(
    tree_file=None,
    nodes_file=None,
    members=None,
    min_support=None,
    min_genomes=2,
    min_length=None,
    prefix_names=None,
    outdir=".",
    prefix="collapsed",
    outgroups=(),
) -> int:
    """Collapse a genome tree into a scheme's hierarchy. Returns an exit status."""
    try:
        tree = read_genome_tree(tree_file, nodes_file)
        # Rooted here as well as in collapse, so that the supports and lengths
        # written next to each node are those of the rooted tree.
        tree = reroot(tree, outgroups)
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
            "removed for its support"
        )

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
        f"wrote {', '.join(written)}"
    )
    return 0


def place_main(
    scheme_dir=None,
    organism=None,
    tree_file=None,
    alignment=None,
    reference=None,
    genome=None,
    members=None,
    min_match=0.5,
    output=None,
    report=None,
) -> int:
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
    logging.info(f"Wrote the genome map to {output}" + (f" and the report to {report}" if report else ""))
    if not placement.clusters:
        logging.error("no genome was placed at any node of the scheme; the genome map is empty")
        return 2
    return 2 if any(finding.severity == PROBLEM for finding in placement.findings) else 0


@dataclass
class CladePlacement:
    """Genomes placed by the clades of the scheme's genomes in a tree."""

    clusters: dict[str, list[str]]
    """The genomes at each node: the scheme's genomes in the tree, by leaf,
    and the genomes placed by the tree."""
    notes: dict[str, tuple[str, str]]
    """The node of every genome and how it was found; an empty node for a
    genome placed nowhere."""


def place_clade(scheme: Scheme, tree_file: str, genome_map: str, outgroups, expected=None) -> CladePlacement:
    """Place the leaves of a tree by the clades of the scheme's genomes.

    Args:
        scheme: The scheme, whose tree file gives the nodes.
        tree_file: The genome tree.
        genome_map: The node of each of the scheme's genomes.
        outgroups: The leaves to root the tree on.
        expected: Names of genomes that should be in the tree, such as the
            records of the alignment it was built from; those the tree lacks
            are reported as placed nowhere rather than left out.

    Raises:
        SchemeError: If the scheme's tree file is not a single rooted tree.
    """
    if scheme.hierarchy is None:
        raise SchemeError("the tree file of the scheme is not a single rooted tree")
    tree = read_genome_tree(tree_file)
    scheme_genomes = {
        genome: node for node, genomes in read_genome_map(genome_map).items() for genome in genomes
    }
    placed, unplaced = place_by_clade(tree, scheme_genomes, scheme.hierarchy.parents, outgroups)

    notes: dict[str, tuple[str, str]] = {}
    clusters: dict[str, list[str]] = {}
    matched, _unmatched, _ambiguous = resolve_names(list(scheme_genomes), tree.leaves)
    outgroup = set(outgroup_leaves(tree, outgroups))
    for genome, leaf in matched.items():
        if leaf not in outgroup:
            notes[leaf] = (scheme_genomes[genome], "scheme genome")
    for leaf, node in placed.items():
        notes[leaf] = (node, "placed by the tree")
    for leaf in unplaced:
        notes[leaf] = ("", "in no clade of the scheme's genomes")
    if expected:
        _in_tree, absent, _ambiguous = resolve_names(list(expected), tree.leaves)
        unplaced_scheme = [name for name in absent if name in scheme_genomes]
        absent = [name for name in absent if name not in scheme_genomes]
        if unplaced_scheme:
            logging.warning(
                f"{len(unplaced_scheme)} scheme genomes of the alignment are not leaves of the "
                f"tree and were placed nowhere: {', '.join(sorted(unplaced_scheme)[:5])}"
            )
        for name in absent:
            notes[name] = ("", "not a leaf of the tree")
        if absent:
            logging.warning(
                f"{len(absent)} genomes are not leaves of the tree and were placed nowhere: "
                f"{', '.join(sorted(absent)[:5])}"
            )
    if unplaced:
        logging.warning(
            f"{len(unplaced)} genomes lie in no clade of the scheme's genomes and were placed "
            f"nowhere: {', '.join(unplaced[:5])}"
        )
    for genome, (node, _note) in notes.items():
        if node:
            clusters.setdefault(node, []).append(genome)
    return CladePlacement(clusters, notes)


def place_clade_main(
    scheme_dir=None,
    organism=None,
    tree_file=None,
    genome_map=None,
    outgroups=(),
    alignment=None,
    output=None,
    report=None,
) -> int:
    """Place genomes by the clades of the scheme's genomes. Returns an exit status."""
    try:
        scheme = read_scheme(scheme_dir, organism)
        expected = fasta_ids(alignment) if alignment else None
        placement = place_clade(scheme, tree_file, genome_map, outgroups, expected)
        with open(output, "w", newline="\n") as handle:
            for node, genomes in sorted(placement.clusters.items()):
                handle.writelines(f"{genome}\t{node}\n" for genome in sorted(genomes))
        if report:
            with open(report, "w", newline="\n") as handle:
                handle.write("genome\tnode\tnote\n")
                handle.writelines(
                    f"{genome}\t{node}\t{note}\n" for genome, (node, note) in sorted(placement.notes.items())
                )
    except (OSError, ValueError, HierarchyError, SchemeError) as error:
        logging.error(error)
        return 1
    logging.info(f"Wrote the genome map to {output}" + (f" and the report to {report}" if report else ""))
    if not placement.clusters:
        logging.error("no genome was placed at any node of the scheme; the genome map is empty")
        return 2
    return 0
