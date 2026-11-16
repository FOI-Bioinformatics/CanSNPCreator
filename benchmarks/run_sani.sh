#!/bin/zsh
# Re-dereplicate an organism's existing genome set at another secondary ANI and
# rerun the scheme and held-out tests.
#   usage: run_sani.sh NAME SRC SANI NODE_PREFIX LABEL REF_GCF OUTGROUP
set -e
NAME=$1; SRC=$2; SANI=$3; PRE=$4; LABEL=$5; REF=$6; OG=$7
B=${0:A:h}; source $B/config.sh
W=$CANSNP_WORK; WD=$W/$NAME
E=$CANSNP_ENVS
REPO=$CANSNP_REPO
export PATH=$E/repgenr_dev/bin:$PATH
PY=$CANSNP_PY
step() { echo "== $NAME: $1 $(date '+%F %T')" }

if [ ! -d $WD/genomes ]; then
  step "ingest from $SRC"
  repgenr ingest -wd $WD --from-workdir $W/$SRC --outgroup $OG
fi
if [ ! -d $WD/derep/representatives ]; then
  step "dereplication at secondary ANI $SANI"
  repgenr dereplicate -wd $WD --tool sourmash --primary-ani 0.9 --secondary-ani $SANI \
    --keeper quality --pin $REF \
  || repgenr dereplicate -wd $WD --tool sourmash --primary-ani 0.9 --secondary-ani $SANI --pin $REF
fi
REFSTEM=$(ls $WD/derep/representatives | grep -v '^\._' | grep "_$REF.fasta" | sed 's/\.fasta$//')
OUTSTEM=$(ls $WD/outgroup | grep -v '^\._' | grep '\.fasta$' | head -1 | sed 's/\.fasta$//')
INFO=$(datasets summary genome accession $REF --as-json-lines 2>/dev/null | $PY -c "
import json, sys
d = json.loads(sys.stdin.readline())
print(d['paired_accession'], d['assembly_info']['assembly_name'])")
GCA=${INFO% *}; ASM=${INFO#* }
P=${REF#GCF_}; URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${P:0:3}/${P:3:3}/${P:6:3}/${REF}_${ASM}/${REF}_${ASM}_genomic.gff.gz
echo "reference $REFSTEM ($GCA, $ASM); outgroup $OUTSTEM; $(ls $WD/derep/representatives | grep -vc '^\._') representatives of $(ls $WD/genomes | grep -vc '^\._')"
step "scheme and held-out test"
$B/bench.sh $NAME $REFSTEM $OUTSTEM $PRE $LABEL $GCA $REF $ASM $URL
step "split-half test"
(cd $REPO && PYTHONPATH=$REPO $PY $B/split.py $W/${NAME}_heldout $REFSTEM $NAME $LABEL $GCA $REF $ASM \
  > $W/${NAME}_heldout/split.log 2>&1) || echo "split.py failed; see ${NAME}_heldout/split.log"
zsh $B/split_summary.sh $W/${NAME}_heldout || true
step finished
