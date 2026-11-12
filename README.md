# CanSNPCreator

CanSNPCreator extracts canonical SNPs (canSNPs) from a whole genome alignment,
scores them for how well they suit a typing assay, and writes the result as the
source files of a [CanSNPer2](https://github.com/FOI-Bioinformatics/CanSNPer2)
scheme.

It also audits a scheme that already exists, reporting where it cannot be typed
with and what would improve it.

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
| alignment | The whole genome alignment in FASTA format, **with the reference genome as its first record**. Every record must be the same length. Either layout is read: gapped, where the reference row carries the aligner's gaps, or projected onto reference coordinates, as x2fa, snippy and RepGenR2 write it, where the reference row has a gap wherever the reference was not aligned. |
| reference | The reference genome in FASTA format, one record per contig, in the order the alignment concatenates them. |
| hierarchy | The scheme as nested YAML or JSON, as a Newick tree, or as a repgenr `child`/`parent` table. |
| group file | The node each genome sits at: a CSV with `Strain` and `Group` columns, a repgenr `representative`/`contained`/`is_self` table, or two tab separated columns of genome and node without a header, as RepGenR2 writes `genomes_map.tsv`. A genome may sit at an internal node. |
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

### Start from a genome tree

A genome tree, such as RepGenR2's `tree/tree.nwk` or `tree2tax.tsv`, has a node
for every split the tree builder found: far more than a scheme has, and many of
them weakly supported. `collapse` removes the internal nodes below a threshold,
joining each into its parent, and places every genome at the deepest node that
remains:

```sh
CanSNPCreator collapse \
    --tree tree.nwk \
    --outgroup Coxiella_sp_GCF_000000000.1 \
    --min-support 0.95 --min-genomes 3 \
    --rename Cb --outdir scheme_draft --prefix coxiella
```

| Option | Removes an internal node when |
| ------ | ----------------------------- |
| `--min-support` | its support is lower, as a fraction; a tree written in percentages is read as fractions |
| `--min-branch-length` | the branch above it is shorter, in the tree's units |
| `--min-genomes` | fewer genomes lie below it (default 2) |

A Newick tree carries its own support values and branch lengths; for
`tree2tax.tsv`, `--nodes-file tree2tax_nodes.tsv` supplies them. Both keep
RepGenR2's outgroup as a leaf unless `tree2tax --remove-outgroup` was given;
its accession is in `outgroup_accession.txt`.

A node that would hold exactly its parent's genomes is joined into the parent
as well, since no SNP could tell the two apart. `--outgroup` names a leaf that
is placed at no node: it stays in the alignment as a genome outside every node,
which is what lets `extract` find the root's canSNP and tell which side of each
split is derived. Without one, `extract` warns that the two sides of the root
cannot be told apart. `--genome-map` reads RepGenR2's `genomes_map.tsv`, so that
genomes folded into a representative by dereplication count towards their
node and are placed with it. `--rename Cb` names the nodes `Cb`, `Cb.1`,
`Cb.1.2` and so on, largest first, in place of the tree's hashes.

It writes `PREFIX_hierarchy.tsv` and `PREFIX_genomes.tsv`, which are the
`--hcfile` and `--groupfile` of `extract`, and `PREFIX_nodes.tsv`, which gives
each node's genome count, support, branch length and name in the tree.

Tree leaves, group files and alignment records are matched by name once
non-alphanumeric characters are removed, and an assembly accession such as
`GCF_000007765.2` matches the one record whose name contains it, since RepGenR2
names records `Family_Genus_species_Accession` but lists accessions alone in
`genomes_map.tsv`.

### Place genomes on an existing scheme

Checking a scheme against genomes needs each genome placed at a node of the
scheme. Typing the genomes with the scheme places them by the SNPs under test,
so a genome that lacks the derived base of its own branch is typed elsewhere
and is never seen to lack it. `place` instead uses each SNP of the scheme only
to identify which clade of a genome tree the node means, and places the
genomes by the tree:

```sh
CanSNPCreator place \
    --scheme-dir CanSNPer2-data/source --organism coxiella_burnetii \
    --tree tree.nwk --fastafile alignment.fasta --reference-file RSA493.fasta \
    -o placed.tsv --report placement.tsv
```

Each node is matched to the side of a branch whose genomes best agree with the
genomes carrying its derived base, by Jaccard index, counting only genomes
with a base at the site. Both sides of every branch are considered, so a tree
rooted at its midpoint, inside a clade of the scheme, is matched as well as one
rooted on an outgroup. A genome is placed at the deepest node whose match
holds it. `--min-match` (default 0.5) sets the lowest index accepted, and the
report gives each node's match, its index and the genomes on either side of
the disagreement. A node with one or two genomes is matched less reliably
than one with many: in a clade of two in which one genome lacks the derived
base, the other genome alone matches as well as the clade.

`audit --tree tree.nwk` places the genomes this way in place of `--groupfile`.

### Audit an existing scheme

```sh
CanSNPCreator audit \
    --scheme-dir path/to/CanSNPer2-data/source \
    --organism francisella_tularensis \
    --reference-file SCHUS4.2.fasta \
    --output findings.tsv
```

Reads a published scheme back in and reports whether it should be extended or
modified. CanSNPer2-data's own `validate.py` answers whether a scheme is well
formed; this answers whether it can be typed with. The checks are:

| Check | What it reports |
| ----- | --------------- |
| tree | a node under two parents, a node with no SNP, a node with two, a SNP naming a node the tree does not hold |
| alleles | a base outside A, C, G and T, which can never be called derived; a derived base equal to the ancestral one; two nodes at one position |
| references | a genome name with an underscore, a malformed accession, a SNP on a genome the references file does not hold |
| callability | nodes with fewer than `--min-required-hits` SNPs on the path from the root, which CanSNPer2 cannot call |
| spacing | two canSNPs closer than `--proximity` bases, where one falls inside the other's assay |
| reference base | with `--reference-file`, a position whose reference base is neither allele, and, when most of them are, a scheme whose coordinates were taken against another version of the assembly |
| better SNPs | with `--score-file` and `--reference-file`, a node where a higher scoring candidate is available, or one the score file does not cover |
| genomes | with `--fastafile`, `--reference-file`, and `--groupfile` or `--tree`, the checks below |

The last group is the only one that tests what defines a canSNP: that every
isolate below the node carries the derived base and no other genome does. Given
the alignment and where its genomes sit in the scheme, the audit reports, for
each SNP of the scheme:

- isolates below the node without the derived base (`not-canonical`), and
  genomes outside it with the derived base (`homoplasy`), by name and by the
  node they sit at;
- SNPs whose isolates carry what the scheme calls ancestral while the genomes
  outside carry what it calls derived (`polarity`);
- genomes with no base at the site (`missing-data`), which is how an accessory
  or poorly aligned region shows;
- for a SNP that fails, other sites that separate exactly that node's isolates
  (`replacement`), or that none does (`unsupported-node`);
- a node that one site alone supports (`single-site`), and a SNP that every
  genome carries, which nothing tests without an outgroup (`no-outgroup`).

It also reads every variable site at which all genomes have one of two bases,
and reports groups of genomes that share SNPs no other genome carries and that
the scheme has no node for (`new-clade`, with the node they would hang from and
how many SNPs mark them), and sites that fit no grouping the scheme can hold
(`conflicting-sites`). A grouping that contradicts a better supported one is
counted as conflicting rather than offered as a clade. Genomes placed at no
node count as outside every node; an outgroup belongs there. On a 2 Mbp
genome, eight genomes take about five seconds; time grows with the number of
genomes.

Each finding is a problem, which makes part of the scheme wrong or impossible
to type, or a recommendation, which leaves it usable. Each also carries the
action it asks for: `modify` an existing SNP, `extend` the scheme with SNPs it
does not have, or `check` something the tool cannot settle itself. A problem
exits 2, a recommendation exits 0 unless `--strict` is given, which is what to
use in continuous integration. `--output` writes the findings as a tab
separated file as well as logging them.

## Exit status

| Status | Meaning |
| ------ | ------- |
| 0 | finished, output written |
| 1 | a file could not be read, or an input is not usable |
| 2 | the inputs are inconsistent with one another, the audited scheme holds a problem, or no subcommand was given; nothing was written |

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

Two tests use the CanSNPer2-data repository when it is checked out beside this
one, and are skipped otherwise: one runs that repository's `validate.py` over
an exported scheme, the other audits its published schemes.

## Contributing

Fork the repository and open a pull request, or open an issue for a bug report
or a feature request.

## License

GPL-3.0-or-later. See `LICENSE`.

## Contact

[Open an issue](https://github.com/FOI-Bioinformatics/CanSNPCreator/issues).
