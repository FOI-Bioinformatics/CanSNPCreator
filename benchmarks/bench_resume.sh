#!/bin/zsh
# Scheme and held-out test of one organism, as bench.sh, but resumed at
# `evaluate`; RepGenR2 metadata, genome and dereplicate must have run in WORK/NAME.
# usage: bench_resume.sh NAME REFERENCE_STEM OUTGROUP_STEM NODE_PREFIX GENOME_LABEL GCA GCF ASSEMBLY GFF_URL
set -e
NAME=$1; REF=$2; OUT=$3; PRE=$4; LABEL=$5; GCA=$6; GCF=$7; ASM=$8; GFFURL=$9
B=${0:A:h}; source $B/config.sh
W=$CANSNP_WORK; WD=$W/$NAME; H=$W/${NAME}_heldout
E=$CANSNP_ENVS
REPO=$CANSNP_REPO
PY=$CANSNP_PY
cs() { (cd $REPO && PATH=$E/mlvatyper/bin:$PATH $PY -m cansnpscorer.cansnp_scorer "$@") }
mkdir -p $H
cs evaluate --snps $H/snps.vcf --reference $H/$REF.fasta --genes $H/reference.gff3 --repeats $H/repeats_trf.bed --hierarchy $H/scheme_hierarchy.tsv -o $H/scores.vcf > $H/evaluate.log 2>&1
rm -rf $H/scheme
cs export --scores $H/scores.vcf --hierarchy $H/scheme_hierarchy.tsv --reference $H/$REF.fasta --outdir $H/scheme --organism ${NAME}_test --genome $LABEL --strain $LABEL --citation NA --genbank-id $GCA --refseq-id $GCF --assembly-name $ASM --genome-map $H/scheme_genomes.tsv --alignment $WD/tree/msa/full_alignment.fasta --snps-per-node 5 --prune > $H/export.log 2>&1

echo "== expected nodes $(date +%T)"
"$PY" - $WD $H $REF <<'PYEOF'
import os, re, sys
WD, H, REF = sys.argv[1:4]
kept = {l.rstrip("\n").split("\t")[-1] for l in open(f"{H}/scheme/" + next(f for f in os.listdir(f"{H}/scheme") if f.endswith("_tree.txt")))}
node = dict(l.rstrip("\n").split("\t")[:2] for l in open(f"{H}/scheme_genomes.tsv"))
def up(x):
    while x not in kept and "." in x:
        x = x.rsplit(".", 1)[0]
    return x
n = moved = 0
with open(f"{H}/heldout_expected.tsv", "w") as out:
    for line in open(f"{WD}/derep/clusters.tsv").readlines()[1:]:
        rep, member = (x.removesuffix(".fasta") for x in line.rstrip("\n").split("\t")[:2])
        if rep == member or rep not in node:
            continue
        moved += up(node[rep]) != node[rep]; n += 1
        out.write(f"{member}\t{up(node[rep])}\n")
    out.write(f"{REF}\t{up(node[REF])}\n")
print(n, "held-out genomes with an expected node;", moved, "moved up from a pruned node")
PYEOF

echo "== audit $(date +%T)"
(/usr/bin/time -l zsh -c "cd $REPO && $PY -m cansnpscorer.cansnp_scorer audit --scheme-dir $H/scheme --reference $H/$REF.fasta --alignment $H/heldout_alignment.fasta --genome-map $H/heldout_expected.tsv --typing-report $H/heldout_typing --report $H/heldout_report.md -o $H/heldout_findings.tsv") > $H/heldout_audit.log 2>&1 || true
echo "== done $(date +%T)"
grep -h -E "Kept|Found [0-9]+ canSNPs|have no canSNP:|Wrote a scheme" $H/collapse.log $H/extract.log $H/export.log | cut -c1-160
cat $H/heldout_typing_summary.tsv | cut -c1-200
