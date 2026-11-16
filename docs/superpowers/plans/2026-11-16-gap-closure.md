# CanSNPCreator gap closure Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the gaps found in the 2026-11-16 SWOT and gap analysis of CanSNPCreator 0.4.0: a reproducible benchmark, a better truth for held-out genomes, lower memory in `audit`, a narrowest-group rule that compares genomes rather than counts, a prune report, and a written rule for multi-SNP calling that CanSNPer2 can adopt.

**Architecture:** Six independent tasks, one branch and one pull request each, in the order below. Tasks 2 and 3 add code to existing modules (`tree_tools.py`, `helpers/file_utils.py`, `typing_report.py`); task 1 adds a `benchmarks/` directory outside the package; task 6 is a document only. No task changes a default or the format of an existing output column.

**Tech Stack:** Python 3.13+, standard library only (no new dependency), pytest, zsh for the benchmark scripts.

**Spec:** the SWOT and gap analysis given in the session of 2026-11-16 (not saved as a file); benchmark results in the memory file `cansnp-ranking-benchmark.md`.

## Global Constraints

- Python 3.13 or later; the package keeps no runtime dependencies.
- `export --snps-per-node` stays at default 1 and `--min-spacing` at 1000 (decision of 2026-11-15).
- No typing subcommand: typing stays in CanSNPer2; CanSNPCreator creates and evaluates schemes (decision of 2026-11-13).
- Existing columns of `_snps.tsv`, `_genomes.tsv` and the VCF outputs keep their names and order; new information goes in new trailing columns.
- Modest scientific language in code, README and commit messages; no Unicode in Nextflow files.
- No Claude attribution lines in commits or pull requests.
- Every pull request passes CI on Python 3.13 and 3.14 before merge.

## Review Focus

- A held-out genome whose name is not a leaf of the tree (absent from the sampled-sites FASTA): expect it reported as unplaced, not silently dropped. Test in Task 2.
- A scheme node none of whose genomes is in the tree: expect the node skipped for placement, and a genome below it placed at the nearest ancestor that is in the tree. Test in Task 2.
- An alignment whose records are not all the same length when read by columns: expect the same `ValueError` text as `read_alignment`. Test in Task 3.
- A node whose canSNPs lack a base in many genomes, so that derived sets differ only by missing data: expect them in one group, as the count rule does. Test in Task 4.
- A `_snps.tsv` written by 0.4.0, without the new column: expect `read_observed` and `evaluate --observed` to work unchanged, using counts. Test in Task 4.

## Corrections to the analysis

- The analysis named a missing warning for short chains under `min_required_hits`. It exists: `check_callability` in `cansnpscorer/cansnp_audit.py:258` reports `min-required-hits` and `unreachable`. No task.
- The analysis named the typing report's per-genome states as the memory cost. They are already one character per canSNP (`GenomeTyping.states`). The cost is `read_alignment`, which keeps every aligned genome whole (S. sonnei: 2.7 GB alignment). Task 3 addresses that.
- `audit --tree` already places genomes on a scheme from a tree, by the Jaccard match of `place`. What the benchmark's `alltree.py` adds is a different rule (deepest scheme node whose representatives' clade holds the genome) and the sampled-sites matrix for building the tree. Task 2 brings both into the package.

---

### Task 1: Benchmark harness in the repository

**Prerequisite:** `/Volumes/sekvens2` mounted (it was unmounted on 2026-11-16 at 19:24 and came back).

**Files:**
- Create: `benchmarks/README.md`, `benchmarks/config.sh`
- Create (copied from `/Volumes/sekvens2/repgenr/work/cansnp_benchmark/`, not the `.bak` files): `benchmarks/bench.sh`, `bench_resume.sh`, `run_org.sh`, `run_sani.sh`, `run_list.sh`, `split_summary.sh`, `split.py`, `project.py`, `compare_published.py`, `mtb_lineage.py`, `tree_derep.py`
- Test: `tests/test_benchmarks.py`

**Interfaces:**
- Produces: `benchmarks/config.sh` exporting `CANSNP_WORK` (default `/Volumes/sekvens2/repgenr/work`), `CANSNP_ENVS` (default `$HOME/miniforge3/envs`), `CANSNP_REPO` (default: the repository root, from the script's own path), `CANSNP_PY` (default `$CANSNP_ENVS/repgenr_dev/bin/python`). Every script sources it instead of its own hard-coded paths.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_benchmarks.py
import pathlib, re, subprocess
BENCH = pathlib.Path(__file__).parent.parent / "benchmarks"

def test_no_hard_coded_paths():
    for path in BENCH.glob("*"):
        if path.suffix in (".sh", ".py"):
            text = path.read_text()
            assert "/Users/" not in text, path
            assert path.name == "config.sh" or "/Volumes/" not in text, path

def test_shell_scripts_parse():
    for path in BENCH.glob("*.sh"):
        subprocess.run(["zsh", "-n", str(path)], check=True)

def test_python_scripts_compile():
    for path in BENCH.glob("*.py"):
        compile(path.read_text(), str(path), "exec")
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `python -m pytest tests/test_benchmarks.py -v`
Expected: FAIL (no `benchmarks/` directory).

- [ ] **Step 3: Copy the scripts, add `config.sh`, replace the hard-coded paths with the four variables**

`alltree.py` is not copied: Task 2 replaces it. `benchmarks/README.md` lists each script with one line of purpose, the order `run_org.sh` -> `bench.sh` -> `split.py` -> `split_summary.sh`, the external tools needed (RepGenR2 in `repgenr_dev`, the `repgenr_snp` env for `--bin-dir simple=`, FastTreeMP, trf from `mlvatyper`, NCBI `datasets`), and the free disk the per-genome alignment scratch takes (1-3 GB per organism).

- [ ] **Step 4: Run the tests**

Run: `python -m pytest tests/test_benchmarks.py -v`
Expected: 3 passed.

- [ ] **Step 5: Run one organism end to end from the repository copy**

Run: `zsh benchmarks/run_sani.sh coxiella_burnetii_check coxiella_burnetii 0.9995 Cb RSA493 GCF_000007765.2 <outgroup stem from coxiella_burnetii/outgroup>`
Expected: `== coxiella_burnetii_check: finished`, and the held-out summary equal to that of `coxiella_burnetii_heldout` (127 held out, single = set 92.1%). Delete `coxiella_burnetii_check*` afterwards.

- [ ] **Step 6: Commit**

```bash
git add benchmarks tests/test_benchmarks.py
git commit -m "Add the benchmark scripts used for the held-out and split-half tests"
```

---

### Task 2: Tree-based expected nodes (`sites` subcommand and `place --rule clade`)

**Files:**
- Create: `cansnpscorer/sites.py`
- Modify: `cansnpscorer/tree_tools.py` (new function after `place`, line 574 onward), `cansnpscorer/cansnp_scorer.py` (new `sites` parser after `place`, line 421 onward; `--rule` on `place` and `audit`)
- Test: `tests/test_sites.py`, `tests/test_tree_tools.py`
- Modify: `README.md`, section "Placing genomes by the tree" (line 529)

**Interfaces:**
- Produces: `sample_sites(snps_vcf: str, reference: Reference, alignments: Sequence[str], sites: int = 20000, seed: int = 0) -> dict[str, str]` in `sites.py`: reads the SNP positions of `snps_vcf`, samples `sites` of them with `random.Random(seed)`, maps each to a column of every alignment with that alignment's own `AlignmentMap`, and returns one row per record; a record in two alignments is taken from the first.
- Produces: `place_by_clade(tree: Hierarchy, scheme_genomes: Mapping[str, str], parents: Mapping[str, str | None], outgroups: Sequence[str]) -> tuple[dict[str, str], list[str]]` in `tree_tools.py`: roots `tree` on the outgroups with `reroot`, finds for each scheme node the clade (the most recent common ancestor) of the genomes at or below it, and returns, for every leaf not in `scheme_genomes`, the deepest scheme node whose clade holds it; second value lists the leaves placed nowhere.
- Produces CLI: `CanSNPCreator sites --snps VCF --reference FASTA --alignment A [--alignment B] [--sites N] [--seed N] -o FASTA`; `place --rule {jaccard,clade}` (default `jaccard`, the current behaviour) with `--genome-map` and `--outgroup` needed for `clade`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_sites.py
def test_two_alignments_with_different_gaps_give_one_column_set(tmp_path):
    # reference ACGTACGT; alignment A has no gaps, alignment B has a gap column
    # after reference base 2; one SNP at reference position 5.
    rows = sample_sites(vcf_with_position(tmp_path, 5), reference, [a_path, b_path], sites=1)
    assert rows["genomeA1"] == "T" and rows["genomeB1"] == "G"
    assert len(set(map(len, rows.values()))) == 1

def test_sampling_is_reproducible(tmp_path):
    first = sample_sites(..., sites=3, seed=7)
    assert first == sample_sites(..., sites=3, seed=7)
```

```python
# tests/test_tree_tools.py
def test_place_by_clade_deepest_node():
    tree = read_genome_tree(newick("((((r1,h1),r2),(r3,h2)),out);"))
    placed, unplaced = place_by_clade(tree, {"r1": "A.1.1", "r2": "A.1", "r3": "A.2"},
                                      {"A": None, "A.1": "A", "A.1.1": "A.1", "A.2": "A"}, ["out"])
    assert placed == {"h1": "A.1.1", "h2": "A.2"} and unplaced == []

def test_place_by_clade_node_without_genomes_in_tree():
    # A.1.1's only genome r1 missing from the tree: h1 goes to A.1
    ...
    assert placed["h1"] == "A.1"

def test_place_by_clade_leaf_outside_every_clade():
    # h3 sits beside the outgroup
    assert "h3" in unplaced
```

- [ ] **Step 2: Run them to make sure they fail**

Run: `python -m pytest tests/test_sites.py tests/test_tree_tools.py -k "sites or clade" -v`
Expected: FAIL (`sample_sites`, `place_by_clade` not defined).

- [ ] **Step 3: Implement `sample_sites` and `place_by_clade`, and the CLI**

`place --rule clade` writes the same genome map as the Jaccard rule (genome, node), so `audit --genome-map` reads it unchanged. `sites` reads each alignment record by record (Task 3's reader when merged, `read_fasta` before), not whole.

- [ ] **Step 4: Run the tests**

Run: `python -m pytest -q`
Expected: all pass.

- [ ] **Step 5: Check against the benchmark**

Run `sites`, FastTreeMP (`-nt -gtr -fastest`), and `place --rule clade` on `escherichia_coli_heldout` with 50,000 sites and seed 0, as `alltree.py` did.
Expected: `placement.tsv` of `alltree.py` and the new genome map give the same node for every held-out genome, or differences only where `alltree.py` and the new code sample different sites; report the count.

- [ ] **Step 6: Commit and open the pull request**

```bash
git commit -m "Place genomes at the deepest scheme node whose clade holds them, and sample SNP sites for the tree"
```

---

### Task 3: Read only the scheme's columns of an alignment in `audit`

**Files:**
- Modify: `cansnpscorer/helpers/file_utils.py` (after `read_alignment`, line 205)
- Modify: `cansnpscorer/cansnp_audit.py:659` (the `read_alignment` call)
- Test: `tests/test_alignment_audit.py`, `tests/test_typing_report.py`

**Interfaces:**
- Produces: `class ColumnRow` with `__getitem__(column: int) -> str` and `__len__() -> int`; reading a column not kept raises `KeyError(f"column {column} was not read")`.
- Produces: `read_alignment_columns(path: str, columns_of: Callable[[str], Collection[int]]) -> dict[str, str | ColumnRow]`: keeps the first record (the reference) whole, calls `columns_of(reference_row)` once to learn the columns, and keeps only those of every later record. Same length and record-count checks and messages as `read_alignment`.
- Consumes: `scheme_columns` (`typing_report.py:289`) to compute the columns from the reference row.

- [ ] **Step 1: Write the failing tests**

```python
def test_column_rows_match_full_read(tmp_path):
    full = read_alignment(path)
    kept = read_alignment_columns(path, lambda reference_row: {0, 4})
    assert kept[first] == full[first]
    assert all(kept[name][4] == full[name][4] for name in full)

def test_unread_column_raises(tmp_path):
    with pytest.raises(KeyError, match="column 2 was not read"):
        read_alignment_columns(path, lambda row: {0})["g1"][2]

def test_ragged_alignment_message_unchanged(tmp_path):
    with pytest.raises(ValueError, match="the records are not aligned"):
        read_alignment_columns(ragged_path, lambda row: {0})

def test_audit_findings_and_typing_unchanged(tmp_path):
    # existing audit fixture with --typing-report, run with both readers
    assert findings_full == findings_columns
    assert genomes_tsv_full == genomes_tsv_columns
```

- [ ] **Step 2: Run them to make sure they fail**

Run: `python -m pytest tests/test_alignment_audit.py tests/test_typing_report.py -k "column" -v`
Expected: FAIL.

- [ ] **Step 3: Implement, and use it in `audit` unless a check needs every column**

`check_against_alignment` (`alignment_audit.py:580`) always ends with `scan_alignment` (line 396), which reads every column to find clades the scheme does not resolve; it runs whenever the genomes are placed (genome map, tree, or the scheme's genomes file). Add `audit --no-new-clades` (default off, so behaviour is unchanged), which skips `scan_alignment` and `check_extension`. `cansnp_audit.main` uses `read_alignment_columns` when `--no-new-clades` is given or when no placement is available (typing report only), and `read_alignment` otherwise; it logs which reader was used.

- [ ] **Step 4: Run the tests**

Run: `python -m pytest -q`
Expected: all pass.

- [ ] **Step 5: Measure**

Run `audit --typing-report --no-new-clades` on `shigella_sonnei_heldout` with `/usr/bin/time -l`, before and after.
Expected: peak resident memory falls well below the alignment size (2.7 GB); identical `_genomes.tsv`.

- [ ] **Step 6: Commit**

```bash
git commit -m "Read only the scheme's columns of an alignment when typing in audit"
```

---

### Task 4: Narrowest-group rule by derived genome sets

**Files:**
- Modify: `cansnpscorer/typing_report.py`: `_snp_rows` (line 455), the `_snps.tsv` header constants (lines 72-85), `Observed` (line 700), `read_observed` (line 716), `width_groups` (line 780), `judge_by_narrowest_group` (line 800)
- Test: `tests/test_typing_report.py`

**Interfaces:**
- Produces: a new trailing column `derived_set` in `_snps.tsv`: per node, canSNPs whose derived genomes are the same set, missing bases allowed to be either state, share an integer label (1 = the set of most genomes); empty where the canSNP has no derived genome.
- Produces: `Observed.derived_set: int = 0` (0 = unknown, as for a 0.4.0 table).
- Changes: `width_groups(rows)` groups by `derived_set` when every row of the node has one, otherwise by counts as now. `judge_by_narrowest_group` is unchanged in signature and logic above `width_groups`.

- [ ] **Step 1: Write the failing tests**

```python
def test_equal_counts_different_genomes_are_two_groups():
    # node with 4 canSNPs: two derived in {g1,g2,g3}, two in {g4,g5,g6}; counts equal
    rows = _snp_rows(...)
    assert labels_of(rows) == [1, 1, 2, 2]

def test_missing_base_does_not_split_a_group():
    # canSNP derived in {g1,g2} with g3 missing, beside two derived in {g1,g2,g3}
    assert labels_of(rows) == [1, 1, 1]

def test_old_table_without_derived_set_reads_and_groups_by_count(tmp_path):
    observed = read_observed(table_without_column)
    assert all(row.derived_set == 0 for row in observed.values())
    assert judge_by_narrowest_group(observed) == expected_from_counts
```

- [ ] **Step 2: Run them to make sure they fail**

Run: `python -m pytest tests/test_typing_report.py -k "derived_set or group" -v`
Expected: FAIL.

- [ ] **Step 3: Implement**

Label the sets in `_snp_rows` from the per-genome state strings already in memory. Two canSNPs share a set when no genome is derived in one and ancestral in the other; assign labels greedily from the canSNP with the most derived genomes. Cost stays linear in genomes per canSNP pair within a node; nodes with more than 200 canSNPs are labelled against the first 200 only, and the rest by count (log the nodes).

- [ ] **Step 4: Run all tests**

Run: `python -m pytest -q`
Expected: all pass.

- [ ] **Step 5: Re-run the split-half test on V. cholerae, M. tuberculosis and Listeria**

Run: `benchmarks/split.py` for the three organisms (Task 1).
Expected: V. cholerae and M. tuberculosis at least 406 of 409 and 552 of 556, as with counts; report Listeria (counts gave 449, the emulated singleton rule 456).

- [ ] **Step 6: Commit**

```bash
git commit -m "Group a node's canSNPs by the genomes they are derived in, not only by their number"
```

---

### Task 5: Report what `--prune` costs

**Files:**
- Modify: `cansnpscorer/cansnp_export.py:419-422` (the prune warning)
- Test: `tests/test_export.py`

**Interfaces:**
- Consumes: the genome map export already reads (`--genome-map`).
- Produces: the prune warning also states how many genomes of the map sat at a pruned node and now sit at an ancestor, and the deepest pruned depth.

- [ ] **Step 1: Write the failing test**

```python
def test_prune_reports_genomes_moved(tmp_path, caplog):
    # three nodes, middle one without canSNP, two genomes mapped to it
    export(..., prune=True, genome_map=map_path)
    assert "2 genomes of the genome map sat at these nodes" in caplog.text
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `python -m pytest tests/test_export.py -k prune_reports -v`
Expected: FAIL.

- [ ] **Step 3: Implement** the count from the genome map; without a genome map the warning is unchanged.

- [ ] **Step 4: Run all tests** — `python -m pytest -q`, all pass.

- [ ] **Step 5: Commit**

```bash
git commit -m "Report how many genomes --prune moves to an ancestor"
```

---

### Task 6: Multi-SNP calling rule for CanSNPer2

**Files:**
- Create: `docs/multi_snp_calling.md`

**Interfaces:**
- Consumes: `set_path`, `set_call`, `node_state` and `derived_off` in `typing_report.py` (lines 207-288), the `_canSNPs.txt` format (README, "Export"), and CanSNPer2 2.2.x's `--min_required_hits` and `--strictness` (PR #53).

- [ ] **Step 1: Write the document**

Sections: the input (`<organism>_canSNPs.txt`, rank 1 = the SNP file's SNP); the state of a node from its canSNPs (majority of called canSNPs, ties and no calls); the call (deepest node on a chain of derived nodes, ties of two chains left uncalled); what `--min_required_hits` and `--strictness` mean under the set rule; the output CanSNPer2 should add per node (derived, ancestral, missing counts); and the benchmark evidence (set rule against single rule on the 20 organisms: genomes typed and agreement, from `cansnp-ranking-benchmark.md`). Each rule cites the function that implements it in CanSNPCreator.

- [ ] **Step 2: Check every rule against the code**

Run: `python -m pytest tests/test_typing_report.py -k set -v` and confirm each rule in the document has a test there; add a test where one is missing.

- [ ] **Step 3: Commit, and open an issue in CanSNPer2 that links the document**

```bash
git add docs/multi_snp_calling.md
git commit -m "Describe the set rule for CanSNPer2 to call nodes from several canSNPs"
```

---

## Order and release

Tasks 1 and 2 first, since Tasks 4 and 5 are checked with the benchmark scripts and Task 2's placement. Task 3 is independent. Task 6 can be written at any time but should follow Task 4, which changes the grouping it describes. Release 0.5.0 after Tasks 2-5 are merged.
