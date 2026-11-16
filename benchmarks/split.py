"""Split-half test of evaluate --observed on held-out genomes.

Half A types every candidate; the candidates are ranked again with what was
observed; five per node are exported from the plain and from the observed
ranking; both schemes type half B.

usage: split.py HELDOUT_DIR REFERENCE_STEM ORGANISM LABEL GCA GCF ASSEMBLY"""
import os
import subprocess
import sys

from cansnpscorer.helpers.file_utils import read_alignment, reference_parser
from cansnpscorer.scheme import read_scheme
from cansnpscorer.typing_report import type_alignment, write_report

D, REFNAME, ORG, LABEL, GCA, GCF, ASM = sys.argv[1:8]
S = f"{D}/split"
os.makedirs(S, exist_ok=True)
REF = f"{D}/{REFNAME}.fasta"
PY = sys.executable
# The fallback file names are the Legionella defaults of the first benchmark.
HIER = f"{D}/scheme_hierarchy.tsv" if os.path.exists(f"{D}/scheme_hierarchy.tsv") else f"{D}/lp_hierarchy.tsv"
GFF = f"{D}/reference.gff3" if os.path.exists(f"{D}/reference.gff3") else f"{D}/GCF_000092545.1.gff3"
META = ["--hierarchy", HIER, "--reference", REF, "--organism", ORG, "--genome", LABEL,
        "--strain", LABEL, "--citation", "NA", "--genbank-id", GCA, "--refseq-id", GCF,
        "--assembly-name", ASM]


def cli(*arguments):
    print("+", " ".join(arguments[:6]), flush=True)
    done = subprocess.run([PY, "-m", "cansnpscorer.cansnp_scorer", *arguments],
                          capture_output=True, text=True)
    if done.returncode not in (0,):
        print(done.stderr[-3000:])
        raise SystemExit(done.returncode)
    return done.stderr


def evaluate(out, observed=None):
    cli("evaluate", "--snps", f"{D}/snps.vcf", "--reference", REF, "--genes",
        GFF, "--repeats",
        f"{D}/repeats_trf.bed", "--hierarchy", HIER, "-o", out,
        *(["--observed", observed] if observed else []))


def export(scores, outdir, *extra):
    cli("export", "--scores", scores, *META, "--outdir", outdir, "--prune", *extra)


if not os.path.exists(f"{S}/scores_plain.vcf"):
    evaluate(f"{S}/scores_plain.vcf")
export(f"{S}/scores_plain.vcf", f"{S}/scheme_all", "--snps-per-node", "all",
       "--min-spacing", "0")
export(f"{S}/scores_plain.vcf", f"{S}/scheme_plain", "--snps-per-node", "5")

alignment = read_alignment(f"{D}/heldout_alignment.fasta")
reference = reference_parser(REF)
names = list(alignment)
first, rest = names[0], sorted(names[1:])
halves = {"A": rest[0::2], "B": rest[1::2]}


def half(which):
    return {first: alignment[first], **{name: alignment[name] for name in halves[which]}}


def typed(scheme_dir, which, prefix):
    report = type_alignment(read_scheme(scheme_dir), half(which), reference)
    write_report(report, prefix)
    return report


typed(f"{S}/scheme_all", "A", f"{S}/all_A")
evaluate(f"{S}/scores_obsA.vcf", f"{S}/all_A_snps.tsv")
export(f"{S}/scores_obsA.vcf", f"{S}/scheme_obsA", "--snps-per-node", "5")
typed(f"{S}/scheme_plain", "B", f"{S}/plain_B")
typed(f"{S}/scheme_obsA", "B", f"{S}/obsA_B")
# The same on the half the observations came from, as the upper bound.
typed(f"{S}/scheme_plain", "A", f"{S}/plain_A")
typed(f"{S}/scheme_obsA", "A", f"{S}/obsA_A")
print("done")
