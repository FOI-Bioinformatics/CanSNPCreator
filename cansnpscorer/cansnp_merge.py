"""The merge subcommand: join schemes built in tiers into one scheme.

A group whose members differ too much for one alignment is built in tiers: a
coarse scheme over the whole group, against one reference, in which one node
holds a subgroup as a leaf, and a finer scheme of that subgroup, against the
same reference, rooted at that node. The node both hold is the seam. Merging
hangs the finer scheme's tree from the seam and joins the two canSNP files.

The node sets are disjoint apart from the seam, so the merge is a union; what
it adds is the checks. Each refusal below is a scheme CanSNPer2 or
CanSNPer2-data would reject later, at a point far from the cause.
"""

import logging
import os
from collections.abc import Sequence
from dataclasses import dataclass

from cansnpscorer.scheme import (
    CANSNP_COLUMNS,
    CANSNP_SUFFIX,
    GENOME_COLUMNS,
    GENOME_ROLES,
    GENOMES_SUFFIX,
    REFERENCE_COLUMNS,
    REFERENCES_SUFFIX,
    TREE_SUFFIX,
    Scheme,
    SchemeError,
    SchemeReference,
    SchemeSnp,
    check_ranks,
    hierarchy_from_paths,
    read_scheme,
)
from cansnpscorer.validate import report

#: A genome listed with two roles keeps the first of these.
ROLE_ORDER = {role: index for index, role in enumerate(GENOME_ROLES)}


@dataclass
class Merged:
    """Schemes joined at their seams, ready to be written."""

    paths: list[list[str]]
    snps: list[SchemeSnp]
    references: list[SchemeReference]
    genomes: list[tuple[str, str, str]]
    seams: list[tuple[str, int]]
    """Each seam, and how many canSNPs of the finer scheme were dropped there."""
    left_out: list[tuple[str, str]]
    """Genome rows left out because their node is no longer a terminal node."""


def reference_key(references: Sequence[SchemeReference]) -> list[tuple[str, ...]]:
    return sorted((r.genome, r.strain, r.genbank_id, r.refseq_id, r.assembly_name)
                  for r in references)


def merge_schemes(schemes: Sequence[Scheme], names: Sequence[str]
                  ) -> tuple[Merged, list[str], list[str]]:
    """Join schemes at their seams. Returns the merge, its problems and its warnings.

    The first scheme is the coarsest. Each later one must share exactly one node
    with those before it, its own root, which is the seam; the seam keeps the
    canSNPs of the scheme before it, the only one that holds a genome outside
    the seam to settle their polarity.
    """
    problems: list[str] = []
    warnings: list[str] = []
    for scheme, name in zip(schemes, names, strict=True):
        problems += [f"{name}: {problem}" for problem in scheme.problems]

    first = schemes[0]
    paths = [list(path) for path in first.paths]
    path_to = {path[-1]: path for path in paths}
    held = set(first.nodes)
    snps = list(first.snp_set)
    seams: list[tuple[str, int]] = []

    for scheme, name in zip(schemes[1:], names[1:], strict=True):
        if scheme.hierarchy is None:
            problems.append(f"{name}: the tree is not a single rooted tree, so it has no seam")
            continue
        root = scheme.hierarchy.root
        shared = set(scheme.nodes) & held
        if not shared:
            problems.append(
                f"{name} shares no node with the schemes before it, so its tree would hang from "
                f"nothing and the merged tree would have two roots; its root {root!r} must be a "
                f"node of the schemes before it")
            continue
        extra = sorted(shared - {root})
        if extra:
            problems.append(
                f"{name} names {len(extra)} nodes that the schemes before it already hold, besides "
                f"its root {root!r}, the seam where the two meet: {', '.join(extra[:5])}. Node "
                f"names must be distinct apart from the seam")
            continue
        for path in scheme.paths:
            if path[0] != root:
                problems.append(f"{name}: the path {';'.join(path)} does not start at the root "
                                f"{root!r}")
                break
            if len(path) > 1:
                joined = path_to[root] + list(path[1:])
                paths.append(joined)
                path_to[joined[-1]] = joined
        dropped = [snp for snp in scheme.snp_set if snp.node == root]
        snps += [snp for snp in scheme.snp_set if snp.node != root]
        seams.append((root, len(dropped)))
        held |= set(scheme.nodes)

    for scheme, name in zip(schemes[1:], names[1:], strict=True):
        if reference_key(scheme.references) != reference_key(first.references):
            problems.append(
                f"{name} names other reference genomes than {names[0]}; the tiers must be built "
                f"against one reference, or their positions do not agree")

    hierarchy, conflicts, tree_problems = hierarchy_from_paths(paths)
    problems += [f"the merged tree: {problem}" for problem in tree_problems]
    for node, parents in conflicts:
        problems.append(f"the merged tree places {node!r} under {' and '.join(parents)}")

    # CanSNPer2 needs every canSNP at a site of its own; tiers extracted
    # independently from one reference can pick one site for two nodes.
    at: dict[tuple[str, int], set[str]] = {}
    for snp in snps:
        at.setdefault((snp.sequence, snp.position), set()).add(snp.node)
    repeated = sorted((site, nodes) for site, nodes in at.items() if len(nodes) > 1)
    if repeated:
        examples = "; ".join(f"{sequence}:{position} by {', '.join(sorted(nodes))}"
                             for (sequence, position), nodes in repeated[:5])
        problems.append(f"{len(repeated)} sites are used by more than one node, and CanSNPer2 "
                        f"needs a site of its own for every canSNP: {examples}")
    check_ranks(snps, "the merged canSNPs", problems)

    genomes, left_out = merge_genomes(schemes, set(hierarchy.leaves) if hierarchy else set(),
                                      problems)
    return (Merged(paths=paths, snps=snps, references=list(first.references), genomes=genomes,
                   seams=seams, left_out=left_out),
            problems, warnings)


def merge_genomes(schemes: Sequence[Scheme], terminal: set[str], problems: list[str]
                  ) -> tuple[list[tuple[str, str, str]], list[tuple[str, str]]]:
    """The genomes of all schemes, at terminal nodes of the merged tree.

    In the coarse scheme the seam is a leaf, so its genomes file lists the finer
    scheme's genomes there; merged, the seam is an inner node and the genomes
    file lists terminal nodes only, so those rows are left out, and the finer
    scheme places the genomes. A genome one scheme uses as its outgroup and
    another places at a node stays placed.
    """
    placed: dict[str, dict[str, list[str]]] = {}
    outgroup: dict[str, str] = {}
    left_out: list[tuple[str, str]] = []
    for scheme in schemes:
        for genome, node, role in scheme.genomes:
            if not node:
                outgroup.setdefault(genome, role)
            elif node in terminal:
                placed.setdefault(genome, {}).setdefault(node, []).append(role)
            else:
                left_out.append((genome, node))

    rows = []
    for genome, nodes in sorted(placed.items()):
        if len(nodes) > 1:
            problems.append(f"{genome} is placed at {len(nodes)} terminal nodes of the merged "
                            f"tree, {', '.join(sorted(nodes))}; a genome sits at one node")
            continue
        (node, roles), = nodes.items()
        rows.append((genome, node, min(roles, key=ROLE_ORDER.__getitem__)))
    moved = {genome for genome, _node in left_out}
    rows += [(genome, "", role) for genome, role in sorted(outgroup.items())
             if genome not in placed and genome not in moved]
    rows.sort(key=lambda row: (ROLE_ORDER[row[2]], row[1], row[0]))
    return rows, sorted(left_out)


def write_merged(merged: Merged, outdir: str, organism: str) -> list[str]:
    """Write the merged scheme's files, as export writes a scheme's."""
    os.makedirs(outdir, exist_ok=True)
    written = []

    tree = os.path.join(outdir, organism + TREE_SUFFIX)
    with open(tree, "w", newline="\n") as handle:
        handle.writelines("\t".join(path) + "\n" for path in merged.paths)
    written.append(tree)

    order = {path[-1]: index for index, path in enumerate(merged.paths)}
    cansnps = os.path.join(outdir, organism + CANSNP_SUFFIX)
    with open(cansnps, "w", newline="\n") as handle:
        handle.write("\t".join(CANSNP_COLUMNS) + "\n")
        handle.writelines("\t".join([
                snp.node, str(snp.rank), snp.reference, snp.genome, snp.sequence,
                str(snp.position), snp.derived, snp.ancestral, snp.score,
                "yes" if snp.typed else "no"]) + "\n" for snp in sorted(merged.snps, key=lambda snp: (order.get(snp.node, len(order)), snp.rank)))
    written.append(cansnps)

    references = os.path.join(outdir, organism + REFERENCES_SUFFIX)
    with open(references, "w", newline="\n") as handle:
        handle.write("\t".join(REFERENCE_COLUMNS) + "\n")
        handle.writelines("\t".join([r.genome, r.strain, r.genbank_id, r.refseq_id,
                                     r.assembly_name]) + "\n" for r in merged.references)
    written.append(references)

    if merged.genomes:
        genomes = os.path.join(outdir, organism + GENOMES_SUFFIX)
        with open(genomes, "w", newline="\n") as handle:
            handle.write("\t".join(GENOME_COLUMNS) + "\n")
            handle.writelines("\t".join(row) + "\n" for row in merged.genomes)
        written.append(genomes)
    return written


def main(schemes: Sequence[str] = (), outdir: str = ".", organism: str | None = None) -> int:
    """Merge schemes given coarsest first. Returns an exit status."""
    if len(schemes) < 2:
        logging.error("merge joins at least two schemes, each given with --scheme")
        return 1
    try:
        read = [read_scheme(directory) for directory in schemes]
    except (OSError, SchemeError, ValueError) as error:
        logging.error(error)
        return 1

    merged, problems, warnings = merge_schemes(read, list(schemes))
    if not report(problems, warnings, "schemes"):
        return 2
    written = write_merged(merged, outdir, organism or read[0].organism)

    for seam, dropped in merged.seams:
        logging.info(f"Joined at the seam {seam!r}: its canSNPs are those of the scheme before "
                     f"it, and the {dropped} of the finer scheme there were dropped")
    if merged.left_out:
        logging.info(
            f"{len(merged.left_out)} genome rows at a node that is no longer terminal were left "
            f"out, since the genomes file lists terminal nodes; the finer scheme places them: "
            f"{', '.join(f'{genome} ({node})' for genome, node in merged.left_out[:5])}")
    nodes = {node for path in merged.paths for node in path}
    logging.info(f"Merged {len(read)} schemes into {len(nodes)} nodes and {len(merged.snps)} "
                 f"canSNPs; wrote {', '.join(written)}")
    return 0
