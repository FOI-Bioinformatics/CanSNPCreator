# CanSNPCreator

CanSNPCreator extracts canonical SNPs (canSNPs) from a whole genome alignment,
scores them for how well they suit a typing assay, and writes the result as the
source files of a [CanSNPer2](https://github.com/FOI-Bioinformatics/CanSNPer2)
scheme.

A canSNP is a substitution carried by every isolate below one node of a
hierarchy and by no other genome, so it marks that branch. Given an alignment,
a hierarchy and the clusters of isolates at its leaves, CanSNPCreator finds
those substitutions, reports which of them sit in well behaved sequence, and
exports a scheme for
[CanSNPer2-data](https://github.com/FOI-Bioinformatics/CanSNPer2-data).

## Installation

Python 3.9 or later.

```sh
pip install git+https://github.com/FOI-Bioinformatics/CanSNPCreator.git
```

The `trf` subcommand also needs Tandem Repeats Finder on the PATH
(`conda install -c bioconda trf`). Nothing else is needed to extract, evaluate
or export.

## Input files

| File | What it holds |
| ---- | ------------- |
| alignment | The whole genome alignment in FASTA format, **with the reference genome as its first record**. Every record must be the same length. |
| reference | The reference genome in FASTA format, one record per contig, in the order the alignment concatenates them. |
| hierarchy | The scheme as nested YAML or JSON, as a Newick tree, or as a repgenr `child`/`parent` table. |
| group file | The isolates of each leaf: either a CSV with `Strain` and `Group` columns, or a repgenr `representative`/`contained`/`is_self` table. |
| GFF | Gene annotation of the reference, used for the gene and mutation terms. |
| TRF | Tandem repeats of the reference, as written by `CanSNPCreator trf`. |

Two things about coordinates are worth knowing before reading any output:

- **Positions are reported by contig, 1-based**, as published schemes record
  them. The alignment's reference row carries the gaps the aligner introduced,
  so a column of the alignment is not a position in the genome; the reference
  file is what the columns are mapped onto.
- **The reference file must be the genome that was aligned.** The tools check
  this base for base and refuse to run otherwise, since a mismatch would move
  every position reported.

Sequence names are matched after non-alphanumeric characters are removed, so
`GCA_000007125.1` in a group file matches `GCA0000071251` in the alignment.
Cluster and node names are matched as written, so a cluster called `A.1` stays
`A.1`.

## Usage

### Extract canSNPs

```sh
CanSNPCreator extract \
    --fastafile alignment.fasta \
    --reference-file reference.fasta \
    --groupfile groups.csv \
    --hcfile hierarchy.yaml \
    -o snps.tsv
```

Writes one row per variable column, with the node it marks and its polarity:

```
contig	position	node	derived_base	ancestral_base
c1	101	root	T	A
c1	201	A	C	G
c2	1	B.2	T	A
c2	101	None	G	A
```

The derived base is the one carried below the node. A row whose node is `None`
is a SNP that marks no branch; it is kept because it still counts as a
neighbour when the others are scored.

### Score them

```sh
CanSNPCreator trf --input reference.fasta --output repeats.trf

CanSNPCreator evaluate \
    --snp-file snps.tsv \
    --reference-file reference.fasta \
    --gff-file genes.gff \
    --trf-file repeats.trf \
    --hcfile hierarchy.yaml \
    --output-file scores.csv
```

Five terms, each 0 to 100 before weighting:

| Term | Meaning |
| ---- | ------- |
| `neighbour_score` | distance to the closest neighbouring SNP; a nearby SNP interferes with an assay |
| `gc_score` | how close the surrounding GC content is to `--gc-goal` |
| `tandem_score` | distance to the closest tandem repeat |
| `gene_score` | whether the SNP falls in a gene, where its context is better defined |
| `mutation_score` | 100 if the substitution is synonymous, 50 if it changes the amino acid, 0 if a stop codon is gained or lost, or if the SNP is outside a gene |

The gene term is counted twice by default, so the highest reachable total is
600. `--score-weight TERM=FACTOR` changes that, for example
`--score-weight gene_score=1`. `--cutoff-value` defaults to 500, which in
practice requires the SNP to fall inside a gene; `--top-n-snp` keeps only the
best few per node.

The run reports which nodes have no canSNP above the cutoff, and which have
fewer than the three CanSNPer2 needs on a path to call a node.

### Export a CanSNPer2 scheme

```sh
CanSNPCreator export \
    --score-file scores.csv \
    --hcfile hierarchy.yaml \
    --reference-file reference.fasta \
    --outdir scheme \
    --organism francisella_tularensis \
    --genome SCHUS4.2 --strain francisella --reference Svensson2009 \
    --genbank-id GCA_000008985.1 --refseq-id GCF_000008985.1 \
    --assembly-name ASM898v1
```

Writes `<organism>_tree.txt`, `<organism>_snp.txt` and
`<organism>_references.txt` in the form CanSNPer2-data expects: one line per
path from the root, exactly one SNP per node (the highest scoring), and
positions counted over the concatenated genome, since the SNP file has no
contig column. The rules of that repository are checked before anything is
written, so the result passes its `validate.py` unchanged.

A node with no canSNP refuses the scheme, because every node of a CanSNPer2
scheme carries one. `--prune` leaves those nodes out instead, re-parenting
their children to the nearest node that has one.

## Exit status

| Status | Meaning |
| ------ | ------- |
| 0 | finished, output written |
| 1 | a file could not be read, or an input is not usable |
| 2 | the inputs are inconsistent with one another, or no subcommand was given; nothing was written |

`--log-level DEBUG` reports more, `--log-level WARNING` only what went wrong.

## Development

```sh
pip install -e ".[dev]"
python -m pytest tests
python -m ruff check cansnpscorer tests
```

The tests build a small dataset with one canSNP planted per node, including a
two-contig reference, a node whose clusters hold two isolates each, a SNP
crossing two branches, a tri-allelic column and a gapped one, and check that
the planted canSNPs are recovered at the planted coordinates with the planted
polarity. CI runs them on Python 3.9, 3.11 and 3.13.

## Contributing

Fork the repository and open a pull request, or open an issue for a bug report
or a feature request.

## License

GPL-3.0-or-later. See `LICENSE`.

## Contact

[Open an issue](https://github.com/FOI-Bioinformatics/CanSNPCreator/issues).
