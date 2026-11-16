#!/bin/zsh
# Scheme and held-out test of one organism, after RepGenR2 metadata, genome and
# dereplicate have run in WORK/NAME.
# usage: bench.sh NAME REFERENCE_STEM OUTGROUP_STEM NODE_PREFIX GENOME_LABEL GCA GCF ASSEMBLY GFF_URL
set -e
NAME=$1; REF=$2; OUT=$3; PRE=$4; LABEL=$5; GCA=$6; GCF=$7; ASM=$8; GFFURL=$9
B=${0:A:h}; source $B/config.sh
W=$CANSNP_WORK; WD=$W/$NAME; H=$W/${NAME}_heldout
E=$CANSNP_ENVS
REPO=$CANSNP_REPO
PY=$CANSNP_PY
cs() { (cd $REPO && PATH=$E/mlvatyper/bin:$PATH $PY -m cansnpscorer.cansnp_scorer "$@") }
mkdir -p $H
ls $WD/derep/representatives | grep -v '^\._' | sed 's/\.fasta$//' > $H/representatives.txt
cp $WD/derep/representatives/$REF.fasta $H/

echo "== alignment and tree $(date +%T)"
if [ ! -s $WD/tree/msa/full_alignment.fasta ] || [ ! -s $WD/tree/tree.nwk -a ! -s $H/tree.nwk ]; then
  PATH=$E/repgenr_dev/bin:$PATH repgenr --bin-dir simple=$E/repgenr_snp/bin phylo -wd $WD --preset species \
    --treebuilder fasttree --msa-source snptype --snptyper simple --reference $REF.fasta -t 10 \
    > $H/repgenr_phylo.log 2>&1 || echo "repgenr phylo did not finish; see $H/repgenr_phylo.log"
fi
if [ -s $WD/tree/tree.nwk ]; then
  cp $WD/tree/tree.nwk $H/tree.nwk; echo "tree: RepGenR2 (FastTree, GTR)" | tee $H/tree_source.txt
elif [ ! -s $H/tree.nwk ]; then
  OMP_NUM_THREADS=10 $E/repgenr_dev/bin/FastTreeMP -nt -boot 1000 -log $H/fasttree_jc.log \
    $WD/tree/msa/core_snp.fasta > $H/fasttree_jc.unrooted.nwk 2> $H/fasttree_jc.err
  $E/repgenr_dev/bin/python - $H/fasttree_jc.unrooted.nwk $H/tree.nwk $OUT <<'PYEOF'
import sys, dendropy
t = dendropy.Tree.get(path=sys.argv[1], schema="newick", preserve_underscores=True)
n = t.find_node_with_taxon_label(sys.argv[3])
t.reroot_at_edge(n.edge, length1=n.edge.length / 2, length2=n.edge.length / 2, update_bipartitions=False)
t.write(path=sys.argv[2], schema="newick", suppress_rooting=True, unquoted_underscores=True)
PYEOF
  echo "tree: FastTreeMP run directly (JC), rooted on $OUT" | tee $H/tree_source.txt
fi

echo "== held-out projection $(date +%T)"
if [ ! -s $H/heldout_alignment.fasta ]; then
  (cd $H && $PY $B/project.py $E/metaquest/bin/minimap2 $REF.fasta $WD/genomes representatives.txt heldout_alignment.fasta 10 > project.log 2>&1)
fi

echo "== scheme $(date +%T)"
curl -s -L "$GFFURL" | gunzip -c > $H/reference.gff3
cs collapse --tree $H/tree.nwk --outgroup $OUT --min-support 0.9 --rename $PRE --outdir $H --prefix scheme > $H/collapse.log 2>&1
cs trf --reference $H/$REF.fasta -o $H/repeats_trf.bed > $H/trf.log 2>&1
(/usr/bin/time -l zsh -c "cd $REPO && $PY -m cansnpscorer.cansnp_scorer extract --alignment $WD/tree/msa/full_alignment.fasta --reference $H/$REF.fasta --hierarchy $H/scheme_hierarchy.tsv --genome-map $H/scheme_genomes.tsv --min-called 0.9 -o $H/snps.vcf") > $H/extract.log 2>&1
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
