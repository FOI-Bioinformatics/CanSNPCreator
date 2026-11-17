# Benchmark scripts

Scripts used for the held-out and split-half tests of the scheme builder.
All paths are set in `config.sh` and can be overridden from the environment:
`CANSNP_WORK` (working directory, default `/Volumes/sekvens2/repgenr/work`),
`CANSNP_ENVS` (conda environments, default `$HOME/miniforge3/envs`),
`CANSNP_REPO` (this repository) and `CANSNP_PY` (Python of the `repgenr_dev`
environment). The scripts are written for zsh.

## Order of use

`run_org.sh` -> `bench.sh` -> `split.py` -> `split_summary.sh`

`run_org.sh` and `run_sani.sh` call the remaining steps in this order.

## Scripts

- `run_org.sh`: select, download and dereplicate one organism with RepGenR2, then run `bench.sh` and `split.py`.
- `run_sani.sh`: re-dereplicate an existing genome set at another secondary ANI and rerun the tests.
- `run_list.sh`: run `run_org.sh` for several organisms in sequence, with one log per organism.
- `bench.sh`: build the scheme and run the held-out test for one organism.
- `bench_resume.sh`: resume `bench.sh` after the tree and scheme steps.
- `split.py`: split-half test of the scheme.
- `split_summary.sh`: summarise the split-half results.
- `project.py`: project genomes onto a reference with minimap2 as one alignment.
- `compare_published.py`: compare calls with those of a published scheme.
- `mtb_lineage.py`: compare calls with the M. tuberculosis complex lineage SNP barcode.
- `tree_derep.py`: choose representatives from a tree over all genomes instead of by ANI.
- `cansnper2_master.sh`: after `bench.sh`, type the held-out genomes with CanSNPer2 itself (`--aligner minimap2`, from a source tree given by `CANSNPER2_SRC`) and with the typing report of this checkout, then compare the calls.
- `compare_cansnper2.py`: compare CanSNPer2 calls with the single-rule, set-rule and expected calls of the typing report.

## External tools

- RepGenR2 in the `repgenr_dev` environment.
- The `repgenr_snp` environment, used for `--bin-dir simple=`.
- FastTreeMP (`repgenr_dev`).
- `trf`, from the `mlvatyper` environment.
- NCBI `datasets`.

## Disk use

The per-genome alignment scratch takes 1 to 3 GB per organism.
