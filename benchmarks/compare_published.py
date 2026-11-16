"""Compare a scheme's calls with those of a published scheme on the same genomes.

Both typing reports must come from the same alignment. For each node of the
new scheme, the published calls of the genomes it calls are counted: a node is
consistent with the published scheme when those calls lie on one line of
descent of the published tree.

usage: compare_published.py new_typing_genomes.tsv published_typing_genomes.tsv published_tree.txt
"""
import re
import sys
from collections import Counter, defaultdict

new_path, old_path, tree_path = sys.argv[1:4]


def calls(path):
    rows = [line.rstrip("\n").split("\t") for line in open(path)]
    at = {name: index for index, name in enumerate(rows[0])}
    return {row[0]: row for row in rows[1:] if row[at["role"]] == "isolate"}, at


new, nat = calls(new_path)
old, oat = calls(old_path)
lineage = {}
for line in open(tree_path):
    path = line.rstrip("\n").split("\t")
    lineage[path[-1]] = path


def related(a, b):
    return a in lineage.get(b, []) or b in lineage.get(a, [])


def major(name):
    found = re.match(r"^(\d\.[A-Z]+)", name)
    return found.group(1) if found else name


both = [g for g in new if new[g][nat["set_call"]] != "NA" and g in old
        and old[g][oat["set_call"]] != "NA"]
print(len(new), "genomes;", len(both), "called by both schemes (set rule);",
      sum(old[g][oat["set_call"]] == "NA" for g in new), "not called by the published one;",
      sum(new[g][nat["set_call"]] == "NA" for g in new), "not called by the new one")
for rule in ("set_call", "single_call"):
    by_node = defaultdict(list)
    for g in both:
        if new[g][nat[rule]] != "NA":
            by_node[new[g][nat[rule]]].append(old[g][oat["set_call"]])
    total = sum(map(len, by_node.values()))
    same = sum(Counter(found).most_common(1)[0][1] for found in by_node.values())
    same_major = sum(Counter(map(major, found)).most_common(1)[0][1]
                     for found in by_node.values())
    on_line = 0
    for found in by_node.values():
        deepest = max(found, key=lambda name: len(lineage.get(name, [])))
        on_line += sum(related(name, deepest) for name in found)
    print(f"new {rule}: {len(by_node)} nodes; genomes sharing the most common published call "
          f"of their node {same}/{total} ({100 * same / total:.1f}%); its main group "
          f"{same_major} ({100 * same_major / total:.1f}%); on one published line of descent "
          f"{on_line} ({100 * on_line / total:.1f}%)")
by_old = defaultdict(Counter)
for g in both:
    by_old[old[g][oat["set_call"]]][new[g][nat["set_call"]]] += 1
print("published call\tgenomes\tnew set calls\tmost common new call (genomes)")
for name, found in sorted(by_old.items(), key=lambda item: -sum(item[1].values()))[:15]:
    call, count = found.most_common(1)[0]
    print(f"{name}\t{sum(found.values())}\t{len(found)}\t{call} ({count})")
