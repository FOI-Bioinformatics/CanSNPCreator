"""Compare a scheme's calls with the published lineage SNP barcode of the
M. tuberculosis complex (TBProfiler's tbdb barcode.bed, H37Rv coordinates).

A lineage is present in a genome when at least half of its barcode SNPs with a
base carry the lineage's allele. The genome's lineage is the deepest present
one whose named ancestors are all present.

For each node of the scheme the lineages of the genomes called at or below it
are counted. The node's lineage is the deepest one that all those genomes share
(a common prefix). A scheme agrees with the barcode when genomes of different
lineages are not mixed below the nodes that are as fine as the barcode.

usage: mtb_lineage.py barcode.bed heldout_alignment.fasta typing_genomes.tsv out.tsv
"""
import sys
from collections import Counter, defaultdict

barcode_path, alignment_path, typing_path, out_path = sys.argv[1:5]

barcode = defaultdict(list)
for line in open(barcode_path):
    fields = line.rstrip("\n").split("\t")
    barcode[fields[3]].append((int(fields[1]), fields[4].upper()))


def parents(name):
    """The named ancestors of a lineage that the barcode holds."""
    parts = name.split(".")
    return [".".join(parts[:i]) for i in range(1, len(parts)) if ".".join(parts[:i]) in barcode]


def lineage(sequence):
    present = set()
    for name, sites in barcode.items():
        called = [sequence[at] for at, _allele in sites if sequence[at] in "ACGT"]
        hits = sum(sequence[at] == allele for at, allele in sites)
        if called and hits * 2 >= len(called):
            present.add(name)
    whole = [name for name in present if all(parent in present for parent in parents(name))]
    if not whole:
        return "NA"
    return max(whole, key=lambda name: (name.count("."), name))


lineages = {}
name, parts = None, []
with open(alignment_path) as handle:
    for line in handle:
        if line.startswith(">"):
            if name is not None:
                lineages[name] = lineage("".join(parts).upper())
            name, parts = line[1:].split()[0], []
        else:
            parts.append(line.strip())
lineages[name] = lineage("".join(parts).upper())

rows = [line.rstrip("\n").split("\t") for line in open(typing_path)]
header, rows = rows[0], rows[1:]
at = {column: header.index(column) for column in header}
with open(out_path, "w") as out:
    out.write("genome\tlineage\tsingle_call\tset_call\n")
    for row in rows:
        out.write(f"{row[0]}\t{lineages.get(row[0], 'NA')}\t{row[at['single_call']]}\t"
                  f"{row[at['set_call']]}\n")

isolates = [row for row in rows if row[at["role"]] == "isolate"]
print(len(isolates), "genomes;", sum(lineages[row[0]] != "NA" for row in isolates),
      "with a barcode lineage;", len({lineages[row[0]] for row in isolates}), "lineages seen")
print("main lineages:", dict(Counter(lineages[row[0]].split(".")[0] for row in isolates)))


def related(a, b):
    return a == b or a.startswith(b + ".") or b.startswith(a + ".")


for rule in ("single_call", "set_call"):
    below = defaultdict(list)
    for row in isolates:
        call = row[at[rule]]
        if call == "NA" or lineages[row[0]] == "NA":
            continue
        parts = call.split(".")
        for depth in range(1, len(parts) + 1):
            below[".".join(parts[:depth])].append(lineages[row[0]])
    mixed = {}
    for node, found in below.items():
        # Genomes below one node should lie on one line of descent of the barcode.
        deepest = max(found, key=lambda name: name.count("."))
        off = sum(not related(name, deepest) for name in found)
        if off and len(set(name.split(".")[0] for name in found)) > 1:
            mixed[node] = (off, len(found))
    # Per genome: is the call a node that holds one main lineage only?
    by_call = defaultdict(Counter)
    for row in isolates:
        if row[at[rule]] != "NA" and lineages[row[0]] != "NA":
            by_call[row[at[rule]]][lineages[row[0]]] += 1
    pure = sum(max(found.values()) for found in by_call.values())
    total = sum(sum(found.values()) for found in by_call.values())
    main_pure = 0
    for found in by_call.values():
        mains = Counter()
        for name, count in found.items():
            mains[name.split(".")[0]] += count
        main_pure += max(mains.values())
    print(f"{rule}: {len(by_call)} nodes called at; genomes sharing the most common lineage "
          f"of their call: {pure} of {total} ({100 * pure / total:.1f}%); sharing its main "
          f"lineage: {main_pure} ({100 * main_pure / total:.1f}%); nodes below which main "
          f"lineages are mixed: {len(mixed)} of {len(below)}")

# The other direction: how finely does the scheme divide each barcode lineage?
by_lineage = defaultdict(Counter)
for row in isolates:
    by_lineage[lineages[row[0]]][row[at["set_call"]]] += 1
print("lineage\tgenomes\tset calls\tmost common set call (genomes)")
for name, found in sorted(by_lineage.items(), key=lambda item: -sum(item[1].values()))[:15]:
    call, count = found.most_common(1)[0]
    print(f"{name}\t{sum(found.values())}\t{len(found)}\t{call} ({count})")
