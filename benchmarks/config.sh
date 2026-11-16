# Shared settings for the benchmark scripts (source this file from zsh).
# Each variable can be overridden from the environment.
_here=${${(%):-%x}:A:h}
export CANSNP_WORK=${CANSNP_WORK:-/Volumes/sekvens2/repgenr/work}
export CANSNP_ENVS=${CANSNP_ENVS:-$HOME/miniforge3/envs}
export CANSNP_REPO=${CANSNP_REPO:-${_here:h}}
export CANSNP_PY=${CANSNP_PY:-$CANSNP_ENVS/repgenr_dev/bin/python}
unset _here
