"""Compare the calls of CanSNPer2 with the typing report of CanSNPCreator.

usage: compare_cansnper2.py CANSNPER2_CALLS TYPING_GENOMES OUT.tsv

CANSNPER2_CALLS has the columns genome and call (NA without a call), as
cansnper2_master.sh writes it. TYPING_GENOMES is the `_genomes.tsv` of
`audit --typing-report`. Only isolates (role `isolate`) in both tables are
compared; the reference and the outgroups are left out, as in the summary.
OUT.tsv lists every genome; the counts are printed.
"""
import csv
import sys
from collections import Counter

calls_path, genomes_path, out_path = sys.argv[1:4]


def read(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


cansnper2 = {row["genome"]: row["call"] or "NA" for row in read(calls_path)}
typing = {row["genome"]: row for row in read(genomes_path) if row["role"] == "isolate"}
genomes = sorted(set(cansnper2) & set(typing))


def relation(call, other):
    """How CanSNPer2's call relates to another call of the same genome."""
    if call == other:
        return "same"
    if call == "NA":
        return "cansnper2-none"
    if other == "NA":
        return "other-none"
    if other.startswith(call + "."):
        return "cansnper2-above"
    if call.startswith(other + "."):
        return "cansnper2-below"
    return "other-branch"


counts = {"single": Counter(), "set": Counter(), "expected": Counter()}
stops = Counter()
with open(out_path, "w", newline="") as handle:
    out = csv.writer(handle, delimiter="\t", lineterminator="\n")
    out.writerow(["genome", "cansnper2", "single_call", "single_stop_reason", "set_call",
                  "expected", "vs_single", "vs_set", "vs_expected"])
    for genome in genomes:
        row = typing[genome]
        call = cansnper2[genome]
        expected = row.get("expected", "NA") or "NA"
        vs = {"single": relation(call, row["single_call"]),
              "set": relation(call, row["set_call"]),
              "expected": relation(call, expected) if expected != "NA" else "NA"}
        for key, value in vs.items():
            counts[key][value] += 1
        if vs["single"] != "same":
            stops[row["single_stop_reason"]] += 1
        out.writerow([genome, call, row["single_call"], row["single_stop_reason"], row["set_call"],
                      expected, vs["single"], vs["set"], vs["expected"]])

n = len(genomes)
typed = sum(cansnper2[genome] != "NA" for genome in genomes)
print(f"genomes compared\t{n}")
print(f"typed by CanSNPer2\t{typed}\t{100 * typed / n:.1f}%" if n else "typed by CanSNPer2\t0")
for key in ("single", "set", "expected"):
    for value, number in sorted(counts[key].items(), key=lambda item: -item[1]):
        print(f"CanSNPer2 vs {key}\t{value}\t{number}\t{100 * number / n:.1f}%")
for reason, number in sorted(stops.items(), key=lambda item: -item[1]):
    print(f"single stop reason where CanSNPer2 differs from single\t{reason or '-'}\t{number}")
