#!/bin/zsh
# Split-half test of evaluate --observed with CanSNPer2 itself typing half B.
#
# usage: split_cansnper2.sh NAME [JOBS]
#
# Runs split.py with this checkout into WORK/NAME_heldout/split_cansnper2,
# then types the isolates of half B with CanSNPer2 (cansnper2_type.sh) under
# the plain scheme and under the scheme ranked with half A's observations,
# and compares the calls with the typing report of split.py
# (compare_cansnper2.py). Writes, per scheme P in plain and obsA:
#   cansnper2_P_B/             CanSNPer2 database, output and calls
#   comparison_P_B.tsv/.txt    CanSNPer2 against the typing report of half B
# The scheme metadata (organism, strain, accessions) is read from the
# held-out scheme.
set -e
NAME=$1; JOBS=${2:-6}
B=${0:A:h}; source $B/config.sh
W=$CANSNP_WORK; WD=$W/$NAME; H=$W/${NAME}_heldout; S=$H/split_cansnper2
REPO=$CANSNP_REPO
PY=${CANSNPER2_PY:-$CANSNP_ENVS/cansnper2/bin/python}
step() { echo "== $NAME: $1 $(date '+%F %T')" }

REFFASTA=$(ls $H/*.fasta | grep -v heldout_alignment | grep -v '/\._' | head -1)
P=$(ls $H/scheme/*_tree.txt | grep -v '/\._'); P=${P%_tree.txt}
IFS=$'\t' read -r LABEL STRAIN GCA GCF ASM < <(sed -n 2p ${P}_references.txt)
mkdir -p $S

step "split.py (CanSNPCreator $(cd $REPO && git describe --tags --always))"
(cd $REPO && SPLIT_DIR=$S PYTHONPATH=$REPO $PY $B/split.py $H ${${REFFASTA:t}%.fasta} $NAME $LABEL $GCA $GCF $ASM) > $S/split.log 2>&1 \
  || { echo "split.py failed, see $S/split.log"; exit 1; }

for p in plain obsA; do
  step "CanSNPer2 on half B, scheme $p"
  awk -F'\t' 'NR>1 && $2=="isolate"{print $1}' $S/${p}_B_genomes.tsv > $S/names_B.txt
  $B/cansnper2_type.sh $S/scheme_$p $REFFASTA $S/names_B.txt $WD/genomes $S/cansnper2_${p}_B $JOBS
  $PY $B/compare_cansnper2.py $S/cansnper2_${p}_B/cansnper2_calls.tsv $S/${p}_B_genomes.tsv \
    $S/comparison_${p}_B.tsv $H/heldout_expected.tsv > $S/comparison_${p}_B.txt
  echo "-- $p"; cat $S/comparison_${p}_B.txt
done
step "done"
