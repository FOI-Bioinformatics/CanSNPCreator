#!/bin/zsh
# The nine organisms of 2026-11-14 evening, one after another, in the order asked.
B=${0:A:h}; source $B/config.sh
run() { $B/run_org.sh "$@" > $CANSNP_WORK/$1_run.log 2>&1; echo "$1 exit $? $(date +%T)" }
run brucella_anthropi Bant ATCC49188 GCF_000017405.1 gtdb Brucella anthropi
run salmonella_paratyphi_a Spa ATCC9150 GCF_000011885.1 taxid:54388 Salmonella enterica GCF_000195995.1 "" Enterobacteriaceae
run treponema_pallidum Tp Nichols GCF_000008605.1 gtdb Treponema pallidum
run shigella_sonnei Ss Ss046 GCF_000092525.1 ncbi Shigella sonnei GCF_000005845.2
run pseudomonas_aeruginosa Pa PAO1 GCF_000006765.1 gtdb Pseudomonas aeruginosa "" 600
run streptococcus_pneumoniae Sp TIGR4 GCF_000006885.1 gtdb Streptococcus pneumoniae "" 600
run klebsiella_pneumoniae Kp HS11286 GCF_000240185.1 gtdb Klebsiella pneumoniae "" 600
run escherichia_coli Ec MG1655 GCF_000005845.2 gtdb Escherichia coli "" 600
run salmonella_typhi St CT18 GCF_000195995.1 taxid:90370 Salmonella enterica GCF_000011885.1 "" Enterobacteriaceae
echo all done
