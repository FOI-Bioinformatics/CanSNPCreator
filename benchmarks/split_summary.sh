#!/bin/zsh
# usage: split_summary.sh HELDOUT_DIR
source ${0:A:h}/config.sh
H=$1; D=$H/split
grep -c "^done" $H/split.log 2>/dev/null || grep -c "^done" $D/split.log
for p in plain_B obsA_B; do
  echo "== $p"
  grep -E "^(genomes|single call is|single call above|single call below|no single call|single rule stopped at a canSNP|single rule stopped at an anc|single rule stopped for too few)" $D/${p}_summary.tsv | cut -f1-3
  awk -F'\t' 'NR==2{print "rank 1 canSNPs disagreeing: "$10}' $D/${p}_ranks.tsv
  awk -F'\t' -v p=$p 'NR==FNR{e[$1]=$2; next} FNR>1 && $2=="isolate"{x=e[$1]; for(k=3;k<=7;k+=4){c=$k; r=(c==x)?"same":(c=="NA"?"none":(index(x,c".")==1?"above":(index(c,x".")==1?"below":"other"))); n[k,r]++}} END{for(k=3;k<=7;k+=4){printf "vs expected, %s:", (k==3?"single":"set"); for(r in n){split(r,a,SUBSEP); if(a[1]==k) printf " %s=%d", a[2], n[r]}; print ""}}' $H/heldout_expected.tsv $D/${p}_genomes.tsv
done
"$CANSNP_PY" - $D <<'PYEOF'
import glob, sys
D = sys.argv[1]
def r1(p):
    return {f[0]: f[3] for f in (l.rstrip("\n").split("\t") for l in open(p)) if f[1] == "1"}
a = r1(glob.glob(D + "/scheme_plain/*_canSNPs.txt")[0]); b = r1(glob.glob(D + "/scheme_obsA/*_canSNPs.txt")[0])
print(len(a), "nodes; rank 1 changed at", sum(a[k] != b.get(k) for k in a))
PYEOF
