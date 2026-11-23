# Changelog

Changes to CanSNPCreator, newest first. Each release has fuller notes, with the
evidence behind them, on its
[GitHub release page](https://github.com/FOI-Bioinformatics/CanSNPCreator/releases).
A version with an **Upgrading** section changes what an earlier command line,
scheme or file does.

## [Unreleased]

### Upgrading

`extract` now stops when no genome of the alignment is placed outside every
node. Such a genome, the outgroup, is what settles which base of a SNP is
ancestral; without one, earlier versions took the reference genome's base as
ancestral and said so only in a warning, so a run given the wrong alignment
still finished and wrote a scheme whose root could be inverted. Name the
outgroup with `--outgroup`, or give `--allow-no-outgroup` to keep the old
behaviour. A run whose genome map already leaves a genome at no node, such as
the reference, goes ahead as before.

### Added

- `extract --outgroup` names the outgroup, and the VCF header records which
  genomes acted as one, as `##outgroup=` (#39).
- `name` names the nodes of a hierarchy after a table of genome labels. A node
  takes a label when every labelled genome below it carries it, and
  `<prefix>_names.tsv` reports how many labelled genomes stand behind each name
  (#40).
- `merge` joins schemes built in tiers into one scheme at the node they share,
  and refuses schemes that disagree about it (#41).

### Fixed

- The conflicting-sites finding of `audit` listed the genomes tied for most
  often involved in an order that changed from run to run, so the same audit
  could name different genomes; ties are now broken by name.

### Development

- CI checks out CanSNPer2 and CanSNPer2-data and runs the eleven tests that
  check exported schemes against them, which had always been skipped there; a
  skipped test now fails the job, and a weekly run catches changes on their
  side (#42).
- The code is formatted with `ruff format`, checked in CI (#43).
- CI fails when line and branch coverage falls below 90% (#44).

## [0.7.1] - 2026-11-22

Includes 0.7.0, which was not released on its own.

### Upgrading

CanSNPCreator writes and reads the canSNP file of a scheme,
`<organism>_canSNPs.txt`, which holds every canSNP of every node with the ones
the scheme types marked, and no longer writes `<organism>_snp.txt`. It is the
format CanSNPer2 2.4.0 reads, and CanSNPer2-data holds; schemes written by this
version need CanSNPer2 2.4.0 or later. `audit` and `place` read a
`<organism>_snp.txt` an earlier version exported as rank 1 of each node, typed,
so an older export is still audited.

### Changed

- `export` always writes the canSNP file, including at `--snps-per-node 1`.
- A genome placed at a node above another kept node is left out of the genomes
  file, with a warning, since CanSNPer2-data's `validate.py` lists the genomes
  of terminal nodes only. It still takes part in finding the canSNPs.
- A genomes file with CanSNPer2-data's `basis` and `source` columns is read.

## [0.6.0] - 2026-11-21

### Upgrading

Version 0.6.0 writes and reads the CanSNPer2 2.3.0 scheme format, in which each
canSNP is a sequence (a GenBank accession) and a position on it. A scheme made
for an earlier version is converted with `CanSNPer2-database --convert`
(CanSNPer2 2.3.0) before `audit` or `place` reads it. A `_snps.tsv` table
written by `audit --typing-report` before 0.6.0 has no `sequence` column and is
refused by `evaluate --observed`; write it again. The `_snps.tsv` of the typing
report has the columns `sequence` and `position` in place of `contig` and
`contig_position`.

### Changed

- `export --contig-aliases` writes RefSeq contig names as their GenBank
  accessions from an NCBI assembly report, and refuses a name without a pair.
- `audit`, the alignment audit and the typing report key canSNPs by sequence
  and position, and map scheme sequence names through `--contig-aliases`.

## [0.5.1] - 2026-11-17

Benchmark scripts only; the package is unchanged from 0.5.0.

### Added

- `benchmarks/cansnper2_master.sh` and `benchmarks/compare_cansnper2.py` type
  held-out genomes with CanSNPer2 and compare its calls with the typing report.

## [0.5.0] - 2026-11-17

### Changed

- The single rule of the typing report follows CanSNPer2 2.2.1 and 2.2.2.
  New option `audit --strictness` (default 0.7), and new stop reasons
  `other-base` and `strictness`.
- `_snps.tsv` gains the trailing columns `derived_set` and `set_apart`; tables
  written by 0.4.0 are still read.

### Added

- `sites` samples SNP columns for building a tree; `place --rule clade` and
  `audit --rule clade` place genomes by clade.
- `audit --no-new-clades` reads only the scheme's columns of the alignment.
- `export --prune` with `--genome-map` reports the genomes moved to an
  ancestor.

## [0.4.0] - 2026-11-15

### Added

- `evaluate` types the genomes and ranks canSNPs by observed behaviour
  (`--observed`); canSNPs in repeated sequence rank last (`--copy-kmer`).
- `audit` types each genome, with `misplaced`, `discordant-genome` and
  `clustered-set` findings and a Markdown report.

### Fixed

- GFF3 features spanning the origin of a circular sequence; name matching, the
  set rule, and several audit, export, extract and tree details.

## [0.3.0] - 2026-11-13

First release under the name CanSNPCreator.

### Changed

- Alignment columns map to 1-based per-contig reference coordinates, so
  positions are right on multi-contig references.
- SNPs and scores are written as VCF 4.3, repeats as BED.
- Python 3.13 or later; no dependencies; one console script, `CanSNPCreator`.

### Added

- `export`, `audit`, `collapse` and `place`, and `extract --min-called`.

[Unreleased]: https://github.com/FOI-Bioinformatics/CanSNPCreator/compare/v0.7.1...HEAD
[0.7.1]: https://github.com/FOI-Bioinformatics/CanSNPCreator/releases/tag/v0.7.1
[0.6.0]: https://github.com/FOI-Bioinformatics/CanSNPCreator/releases/tag/v0.6.0
[0.5.1]: https://github.com/FOI-Bioinformatics/CanSNPCreator/releases/tag/v0.5.1
[0.5.0]: https://github.com/FOI-Bioinformatics/CanSNPCreator/releases/tag/v0.5.0
[0.4.0]: https://github.com/FOI-Bioinformatics/CanSNPCreator/releases/tag/v0.4.0
[0.3.0]: https://github.com/FOI-Bioinformatics/CanSNPCreator/releases/tag/v0.3.0
