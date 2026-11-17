#!/bin/zsh
# Type the held-out genomes of one organism with CanSNPer2 itself, and again
# with the typing report of this CanSNPCreator checkout, then compare the two.
#
# usage: cansnper2_master.sh NAME [JOBS]
#
# Needs a finished bench.sh run (WORK/NAME and WORK/NAME_heldout). Writes
# WORK/NAME_heldout/cansnper2_master/:
#   NAME.db, refs/        CanSNPer2 database built from the exported scheme
#   calls_N/              CanSNPer2 output, one folder per parallel job
#   cansnper2_calls.tsv   genome and CanSNPer2 call (NA without one)
#   typing_*              the typing report of this checkout (audit, genome map)
#   comparison.tsv/.txt   compare_cansnper2.py
# The held-out alignment is projected again when it is missing.
#
# CanSNPer2 is run from CANSNPER2_SRC (a copy of its source tree, for example
# from `git archive`) with --aligner minimap2 and its default
# --min_required_hits 3 and --strictness 0.7. CANSNPER2_PY must be Python
# 3.13 or later, since it also runs CanSNPCreator.
set -e
NAME=$1; JOBS=${2:-6}
B=${0:A:h}; source $B/config.sh
W=$CANSNP_WORK; WD=$W/$NAME; H=$W/${NAME}_heldout; C=$H/cansnper2_master
E=$CANSNP_ENVS; REPO=$CANSNP_REPO
SRC=${CANSNPER2_SRC:-$W/cansnp_benchmark/cansnper2_master}
PY=${CANSNPER2_PY:-$E/cansnper2/bin/python}
MM=${CANSNP_MINIMAP2:-$E/repgenr_dev/bin/minimap2}
step() { echo "== $NAME: $1 $(date '+%F %T')" }

REFFASTA=$(ls $H/*.fasta | grep -v heldout_alignment | head -1)
SCHEME=$(ls $H/scheme/*_tree.txt); P=${SCHEME%_tree.txt}
REFGENOME=$(awk -F'\t' 'NR==2{print $1}' ${P}_references.txt)
mkdir -p $C/refs $C/logs

if [[ ! -s $H/heldout_alignment.fasta ]]; then
  step "project"
  (cd $H && $PY $B/project.py $MM ${REFFASTA:t} $WD/genomes representatives.txt heldout_alignment.fasta 10 > project.log 2>&1)
fi

step "audit (CanSNPCreator $(cd $REPO && git describe --tags --always))"
(cd $REPO && /usr/bin/time -l $PY -m cansnpscorer.cansnp_scorer audit --scheme-dir $H/scheme \
   --reference $REFFASTA --alignment $H/heldout_alignment.fasta --genome-map $H/heldout_expected.tsv \
   --typing-report $C/typing --report $C/report.md -o $C/findings.tsv) > $C/audit.log 2>&1 || true
# audit exits non-zero when it reports findings; the typing report must exist
[[ -s $C/typing_genomes.tsv ]] || { echo "audit failed, see $C/audit.log"; exit 1; }

step "database"
ln -sf $REFFASTA $C/refs/$REFGENOME.fna
rm -f $C/$NAME.db
(cd $C && PYTHONPATH=$SRC $PY -m CanSNPer2.SNPDatabase -db $C/$NAME.db --tree ${P}_tree.txt \
   --annotation ${P}_snp.txt --references ${P}_references.txt --create) > $C/database.log 2>&1

step "CanSNPer2 on $(wc -l < $H/heldout_expected.tsv | tr -d ' ') genomes, $JOBS jobs"
rm -rf $C/calls_*(N) $C/tmp_*(N)
cut -f1 $H/heldout_expected.tsv | sed "s|^|$WD/genomes/|;s|\$|.fasta|" > $C/queries.txt
split -n r/$JOBS $C/queries.txt $C/queries_ 2>/dev/null || {
  # BSD split has no -n: deal the queries round robin
  awk -v n=$JOBS -v p=$C/queries_ '{print > (p sprintf("%02d", (NR-1)%n))}' $C/queries.txt
}
for part in $C/queries_*; do
  k=${part:t}
  (PYTHONPATH=$SRC $PY -m CanSNPer2.CanSNPerTree -db $C/$NAME.db --refdir $C/refs \
     --aligner minimap2 --aligner_path $MM -o $C/calls_$k --tmpdir $C/tmp_$k --logdir $C/logs \
     --summary --keep_going $(cat $part) > $C/logs/$k.log 2>&1; rm -rf $C/tmp_$k) &
done
wait
rm -f $C/queries_*(N)
cat $C/calls_*/snps_summary.txt | awk -F'\t' 'BEGIN{OFS="\t"; print "genome","call"} {print $1, ($2==""?"NA":$2)}' > $C/cansnper2_calls.tsv

step "compare"
$PY $B/compare_cansnper2.py $C/cansnper2_calls.tsv $C/typing_genomes.tsv $C/comparison.tsv | tee $C/comparison.txt
step "done"
