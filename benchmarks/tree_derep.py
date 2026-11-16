"""Choose representatives from a tree over every genome, not by ANI.

Where every genome of a species lies above the resolution of an ANI tool, the
clusters are wider than the clades and the scheme loses the groups a published
scheme resolves. This cuts a tree over all genomes at the height that gives the
wanted number of clusters and writes the result as a RepGenR2 derep directory.

usage: tree_derep.py TREE WORKDIR OUTGROUP_STEM TARGET [PIN_ACCESSION]
"""
import json
import os
import shutil
import sys

import dendropy

TREE, WD, OUTGROUP, TARGET = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
PIN = sys.argv[5] if len(sys.argv) > 5 else None

tree = dendropy.Tree.get(path=TREE, schema="newick", preserve_underscores=True)
node = tree.find_node_with_taxon_label(OUTGROUP)
if node is None:
    raise SystemExit(f"the outgroup {OUTGROUP} is not in the tree")
tree.reroot_at_edge(node.edge, update_bipartitions=True)
tree.prune_taxa_with_labels([OUTGROUP])
tree.calc_node_ages(ultrametricity_precision=float("inf"))

genomes = {name[:-6] for name in os.listdir(f"{WD}/genomes")
           if name.endswith(".fasta") and not name.startswith("._")}
tips = [leaf for leaf in tree.leaf_node_iter() if leaf.taxon.label in genomes]
print(f"{len(tips)} of the tree's tips are genomes of the working directory")


distances = tree.phylogenetic_distance_matrix()
taxa = [leaf.taxon for leaf in tips]


def clusters(limit):
    """Tips grouped so that each lies within limit of its group's first tip."""
    heads, groups = [], []
    for taxon in taxa:
        near = None
        for index, head in enumerate(heads):
            if distances.distance(taxon, head) <= limit:
                near = index
                break
        if near is None:
            heads.append(taxon)
            groups.append([taxon.label])
        else:
            groups[near].append(taxon.label)
    return groups


low, high = 0.0, max(distances.distance(taxa[0], taxon) for taxon in taxa[1:])
for _ in range(18):
    middle = (low + high) / 2
    if len(clusters(middle)) > TARGET:
        low = middle
    else:
        high = middle
groups = clusters(high)
print(f"{len(groups)} clusters at a distance of {high:.3g}")

# One representative per cluster: the pinned genome where it is in the cluster,
# otherwise the longest assembly.
sizes = {}
for labels in groups:
    for label in labels:
        path = f"{WD}/genomes/{label}.fasta"
        sizes[label] = os.path.getsize(os.path.realpath(path))
chosen = []
for labels in groups:
    pinned = [l for l in labels if PIN and PIN in l]
    chosen.append((pinned[0] if pinned else max(labels, key=lambda l: sizes[l]), labels))

derep = f"{WD}/derep"
if os.path.exists(derep):
    shutil.move(derep, f"{derep}_ani")
os.makedirs(f"{derep}/representatives")
with open(f"{derep}/clusters.tsv", "w") as handle:
    handle.write("representative\tmember\n")
    for representative, labels in chosen:
        os.symlink(os.path.realpath(f"{WD}/genomes/{representative}.fasta"),
                   f"{derep}/representatives/{representative}.fasta")
        for label in labels:
            handle.write(f"{representative}.fasta\t{label}.fasta\n")
with open(f"{derep}/derep_stats.json", "w") as handle:
    json.dump({"tool": "tree_derep.py", "tree": TREE, "height": high,
               "representatives": len(chosen), "genomes": len(sizes)}, handle, indent=2)
print(f"wrote {len(chosen)} representatives to {derep}/representatives")
