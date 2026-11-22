#!/bin/zsh
# Type genomes with CanSNPer2 itself, from a scheme exported by CanSNPCreator.
#
# usage: cansnper2_type.sh SCHEME_DIR REFFASTA NAMES GENOME_DIR OUT [JOBS]
#
# SCHEME_DIR holds the exported *_tree.txt, *_canSNPs.txt and *_references.txt.
# NAMES lists one genome per line; GENOME_DIR holds NAME.fasta or
# PREFIX_NAME.fasta. Each query is linked into OUT/genomes under its accession
# (GCA_/GCF_ number) when the name holds one, so that CanSNPer2 names the
# sample as compare_cansnper2.py expects. Writes OUT/cansnper2_calls.tsv
# (genome and call, NA without one), with CanSNPer2 output in OUT/calls_*.
#
# CanSNPer2 is run from CANSNPER2_SRC with --aligner minimap2 and its default
# --min_required_hits 3 and --strictness 0.7.
set -e
SCHEME=$1; REFFASTA=$2; NAMES=$3; GDIR=$4; OUT=$5; JOBS=${6:-6}
B=${0:A:h}; source $B/config.sh
W=$CANSNP_WORK; E=$CANSNP_ENVS
SRC=${CANSNPER2_SRC:-$W/cansnp_benchmark/cansnper2_master}
PY=${CANSNPER2_PY:-$E/cansnper2/bin/python}
MM=${CANSNP_MINIMAP2:-$E/repgenr_dev/bin/minimap2}

P=$(ls $SCHEME/*_tree.txt | grep -v '/\._'); P=${P%_tree.txt}
REFGENOME=$(awk -F'\t' 'NR==2{print $1}' ${P}_references.txt)
rm -rf $OUT/refs $OUT/genomes $OUT/calls_*(N) $OUT/tmp_*(N) $OUT/queries_*(N)
mkdir -p $OUT/refs $OUT/genomes $OUT/logs
ln -s $REFFASTA $OUT/refs/$REFGENOME.fna

rm -f $OUT/scheme.db
(cd $OUT && PYTHONPATH=$SRC $PY -m CanSNPer2.SNPDatabase -db $OUT/scheme.db --tree ${P}_tree.txt \
   --annotation ${P}_canSNPs.txt --references ${P}_references.txt --create) > $OUT/database.log 2>&1

while read -r g; do
  [[ -n $g ]] || continue
  f=$GDIR/$g.fasta; [[ -e $f ]] || f=$(ls $GDIR/*_$g.fasta 2>/dev/null | head -1)
  [[ -n $f && -e $f ]] || { echo "no genome file for $g in $GDIR"; exit 1; }
  short=$(echo $g | grep -oE 'GC[AF]_[0-9]+\.[0-9]+$' || echo $g)
  ln -s $f $OUT/genomes/$short.fasta
done < $NAMES
ls $OUT/genomes/*.fasta > $OUT/queries.txt
# deal the queries round robin (BSD split has no -n)
awk -v n=$JOBS -v p=$OUT/queries_ '{print > (p sprintf("%02d", (NR-1)%n))}' $OUT/queries.txt
for part in $OUT/queries_*; do
  k=${part:t}
  (PYTHONPATH=$SRC $PY -m CanSNPer2.CanSNPerTree -db $OUT/scheme.db --refdir $OUT/refs \
     --aligner minimap2 --aligner_path $MM -o $OUT/calls_$k --tmpdir $OUT/tmp_$k --logdir $OUT/logs \
     --summary --keep_going $(cat $part) > $OUT/logs/$k.log 2>&1; rm -rf $OUT/tmp_$k) &
done
wait
rm -f $OUT/queries_*(N)
# snps_summary.txt from CanSNPer2 2.3.0 starts with # lines on the run and a header row
# (sample final_snp support snp_path); earlier versions wrote the rows only. Both are read.
cat $OUT/calls_*/snps_summary.txt | awk -F'\t' 'BEGIN{OFS="\t"; print "genome","call"}
  /^#/ || ($1=="sample" && $2=="final_snp") {next}
  {print $1, ($2==""?"NA":$2)}' > $OUT/cansnper2_calls.tsv
