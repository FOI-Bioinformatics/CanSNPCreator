"""Project genomes onto a reference with minimap2, as one alignment in
reference coordinates: the reference first, then one record per genome, each
as long as the reference. A base is the genome's where exactly one primary
alignment of mapping quality 20 or more covers the position, '-' where none
does or the genome has a deletion, and 'N' where more than one does.

usage: project.py minimap2 reference.fasta genomes_dir skip_list out.fasta [threads]
"""
import os
import re
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor

minimap2, reference_path, genomes_dir, skip_path, out_path = sys.argv[1:6]
threads = int(sys.argv[6]) if len(sys.argv) > 6 else 8
CS = re.compile(r":(\d+)|\*([a-z])([a-z])|\+([a-z]+)|-([a-z]+)")


def read_fasta(path):
    records, name, parts = [], None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(parts).upper()))
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    records.append((name, "".join(parts).upper()))
    return records


contigs = read_fasta(reference_path)
offsets, total = {}, 0
for name, sequence in contigs:
    offsets[name] = total
    total += len(sequence)
reference = "".join(sequence for _name, sequence in contigs)


def project(path):
    paf = subprocess.run([minimap2, "-cx", "asm20", "--cs", "-t", "1", reference_path, path],
                         capture_output=True, text=True, check=True).stdout
    out = bytearray(b"-" * total)
    depth = bytearray(total)
    for line in paf.splitlines():
        fields = line.split("\t")
        if int(fields[11]) < 20 or "tp:A:P" not in fields[12:]:
            continue
        cs = next(tag[5:] for tag in fields[12:] if tag.startswith("cs:Z:"))
        at = offsets[fields[5]] + int(fields[7])
        for match, _ref, query, _insertion, deletion in CS.findall(cs):
            if match:
                n = int(match)
                out[at:at + n] = reference[at:at + n].encode()
                for i in range(at, at + n):
                    depth[i] += depth[i] < 2
                at += n
            elif query:
                out[at] = ord(query.upper())
                depth[at] += depth[at] < 2
                at += 1
            elif deletion:
                for i in range(at, at + len(deletion)):
                    depth[i] += depth[i] < 2
                at += len(deletion)
    for i in range(total):
        if depth[i] > 1:
            out[i] = 78  # N
    return os.path.basename(path).removesuffix(".fasta"), out.decode()


def main():
    with open(skip_path) as handle:
        skip = {line.split()[0] for line in handle if line.strip()}
    paths = sorted(os.path.join(genomes_dir, name) for name in os.listdir(genomes_dir)
                   if name.endswith(".fasta") and not name.startswith("._")
                   and name.removesuffix(".fasta") not in skip)
    print(f"{len(paths)} genomes to project onto {total} bases", flush=True)
    with open(out_path, "w") as handle, ProcessPoolExecutor(threads) as pool:
        handle.write(f">{os.path.basename(reference_path).removesuffix('.fasta')}\n{reference}\n")
        for done, (name, sequence) in enumerate(pool.map(project, paths, chunksize=4), start=1):
            handle.write(f">{name}\n{sequence}\n")
            if done % 50 == 0:
                print(f"{done} of {len(paths)}", flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
