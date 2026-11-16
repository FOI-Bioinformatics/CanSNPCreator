#!/bin/zsh
# Select, download and dereplicate one organism with RepGenR2, then run
# bench.sh (scheme and held-out test) and split.py (split-half test).
#
# usage: run_org.sh NAME NODE_PREFIX LABEL REF_GCF MODE GENUS SPECIES [OUTGROUP_GCF] [LIMIT] [FAMILY]
#   MODE gtdb      GTDB r232 species (outgroup: the sister representative)
#   MODE ncbi      NCBI species (OUTGROUP_GCF required)
#   MODE taxid:N   every current RefSeq assembly of NCBI taxon N (any rank),
#                  downloaded and ingested; OUTGROUP_GCF and FAMILY required
# LABEL must not contain an underscore.
set -e
NAME=$1; PRE=$2; LABEL=$3; REF=$4; MODE=$5; G=$6; S=$7; OG=${8:-}; LIMIT=${9:-}; FAM=${10:-}
B=${0:A:h}; source $B/config.sh
W=$CANSNP_WORK; WD=$W/$NAME
E=$CANSNP_ENVS
REPO=$CANSNP_REPO
export PATH=$E/repgenr_dev/bin:$PATH
PY=$CANSNP_PY
step() { echo "== $NAME: $1 $(date '+%F %T')" }

# One assembly as DIR/STEM_ACC.fasta.
fetch() {
  local acc=$1 dir=$2 stem=$3 tmp=$W/.fetch_$1
  rm -rf $tmp; mkdir -p $tmp $dir
  datasets download genome accession $acc --include genome --filename $tmp/d.zip > /dev/null
  (cd $tmp && unzip -q d.zip)
  cat $tmp/ncbi_dataset/data/$acc/*.fna > $dir/${stem}_$acc.fasta
  rm -rf $tmp
}

if [ ! -d $WD/genomes ]; then
  step selection
  LIM=(); [ -n "$LIMIT" ] && LIM=(--limit $LIMIT)
  case $MODE in
    gtdb)
      repgenr metadata -wd $WD -d all -l species -r 232.0 --gtdb-version bac120 \
        --metadata-path $W/../live-cache/bac120_metadata_r232.tsv.gz -tg $G -ts $S $LIM
      repgenr genome -wd $WD ;;
    ncbi)
      repgenr metadata -wd $WD -d all -l species --source ncbi -tg $G -ts $S \
        --outgroup-accession $OG $LIM
      repgenr genome -wd $WD ;;
    taxid:*)
      DL=$W/${NAME}_download; mkdir -p $DL/genomes $DL/outgroup
      if [ ! -s $DL/d.zip ]; then
        datasets download genome taxon ${MODE#taxid:} --assembly-source refseq \
          --assembly-version latest --exclude-atypical --include genome --filename $DL/d.zip
      fi
      (cd $DL && rm -rf ncbi_dataset && unzip -q -o d.zip)
      for d in $DL/ncbi_dataset/data/GC[AF]_*(/); do
        cat $d/*.fna > $DL/genomes/${FAM}_${G}_${S}_${d:t}.fasta
      done
      rm -rf $DL/ncbi_dataset
      [ -s $DL/outgroup/${FAM}_${G}_outgroup_$OG.fasta ] || fetch $OG $DL/outgroup ${FAM}_${G}_outgroup
      repgenr ingest -wd $WD --genomes-dir $DL/genomes --outgroup $DL/outgroup/${FAM}_${G}_outgroup_$OG.fasta ;;
  esac
fi

# --pin needs the reference in the workdir; a capped selection may leave it out.
if ! ls $WD/genomes | grep -v '^\._' | grep -q "_$REF.fasta"; then
  step "adding the reference $REF"
  FIRST=$(ls $WD/genomes | grep -v '^\._' | grep fasta | head -1)
  fetch $REF $W/${NAME}_reference ${FIRST%_GC[AF]_*}
  OGF=$(ls $WD/outgroup/*.fasta | grep -v '/\._' | head -1)
  mv $WD ${WD}_selection
  repgenr ingest -wd $WD --from-workdir ${WD}_selection --genomes-dir $W/${NAME}_reference \
    --outgroup ${OGF/$WD/${WD}_selection}
fi

if [ ! -d $WD/derep/representatives ]; then
  step dereplication
  repgenr dereplicate -wd $WD --tool sourmash --primary-ani 0.9 --secondary-ani 0.9995 \
    --keeper quality --pin $REF \
  || repgenr dereplicate -wd $WD --tool sourmash --primary-ani 0.9 --secondary-ani 0.9995 --pin $REF
fi

REFSTEM=$(ls $WD/derep/representatives | grep -v '^\._' | grep "_$REF.fasta" | sed 's/\.fasta$//')
OUTSTEM=$(ls $WD/outgroup | grep -v '^\._' | grep '\.fasta$' | head -1 | sed 's/\.fasta$//')
INFO=$(datasets summary genome accession $REF --as-json-lines 2>/dev/null | $PY -c "
import json, sys
d = json.loads(sys.stdin.readline())
print(d['paired_accession'], d['assembly_info']['assembly_name'])")
GCA=${INFO% *}; ASM=${INFO#* }
P=${REF#GCF_}; URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${P:0:3}/${P:3:3}/${P:6:3}/${REF}_${ASM}/${REF}_${ASM}_genomic.gff.gz
echo "reference $REFSTEM ($GCA, $ASM); outgroup $OUTSTEM; $(ls $WD/derep/representatives | grep -vc '^\._') representatives"

step "scheme and held-out test"
$B/bench.sh $NAME $REFSTEM $OUTSTEM $PRE $LABEL $GCA $REF $ASM $URL
step "split-half test"
(cd $REPO && PYTHONPATH=$REPO $PY $B/split.py $W/${NAME}_heldout $REFSTEM $NAME $LABEL $GCA $REF $ASM \
  > $W/${NAME}_heldout/split.log 2>&1) || echo "split.py failed; see ${NAME}_heldout/split.log"
zsh $B/split_summary.sh $W/${NAME}_heldout || true
step finished
