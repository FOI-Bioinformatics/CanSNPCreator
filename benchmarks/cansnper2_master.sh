#!/bin/zsh
# Type the held-out genomes of one organism with CanSNPer2 itself, and again
# with the typing report of this CanSNPCreator checkout, then compare the two.
#
# usage: cansnper2_master.sh NAME [JOBS]
#
# Needs a finished bench.sh run (WORK/NAME and WORK/NAME_heldout). Writes
# WORK/NAME_heldout/cansnper2_master/:
#   cansnper2_calls.tsv   genome and CanSNPer2 call (NA without one), from
#                         cansnper2_type.sh, with its database and output
#   typing_*              the typing report of this checkout (audit, genome map)
#   comparison.tsv/.txt   compare_cansnper2.py
# The held-out alignment is projected again when it is missing.
#
# The held-out scheme is read as the CanSNPer2 2.3.0 format, so CANSNPER2_SRC
# must be a 2.3.0 source tree and cansnper2_type.sh must create the database
# with 2.3.0.
# CanSNPer2 is run from CANSNPER2_SRC (a copy of its source tree, for example
# from `git archive`) with --aligner minimap2 and its default
# --min_required_hits 3 and --strictness 0.7. CANSNPER2_PY must be Python
# 3.13 or later, since it also runs CanSNPCreator.
set -e
NAME=$1; JOBS=${2:-6}
B=${0:A:h}; source $B/config.sh
W=$CANSNP_WORK; WD=$W/$NAME; H=$W/${NAME}_heldout; C=$H/cansnper2_master
E=$CANSNP_ENVS; REPO=$CANSNP_REPO
PY=${CANSNPER2_PY:-$E/cansnper2/bin/python}
MM=${CANSNP_MINIMAP2:-$E/repgenr_dev/bin/minimap2}
step() { echo "== $NAME: $1 $(date '+%F %T')" }

REFFASTA=$(ls $H/*.fasta | grep -v heldout_alignment | head -1)
mkdir -p $C

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

step "CanSNPer2 on $(wc -l < $H/heldout_expected.tsv | tr -d ' ') genomes, $JOBS jobs"
cut -f1 $H/heldout_expected.tsv > $C/names.txt
$B/cansnper2_type.sh $H/scheme $REFFASTA $C/names.txt $WD/genomes $C $JOBS

step "compare"
$PY $B/compare_cansnper2.py $C/cansnper2_calls.tsv $C/typing_genomes.tsv $C/comparison.tsv | tee $C/comparison.txt
step "done"
