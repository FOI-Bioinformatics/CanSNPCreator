# CanSNPCreator

CanSNPCreator finds canonical SNPs (canSNPs) in a whole genome alignment,
scores them for how well they suit a typing assay, and writes the result as the
source files of a [CanSNPer2](https://github.com/FOI-Bioinformatics/CanSNPer2)
scheme. It also audits a scheme that already exists, reporting where it cannot
be typed with and how it would be extended or modified.

A canSNP is a substitution carried by every genome below one node of a scheme
and by no other genome, so it marks that branch. Given an alignment, the
scheme's tree and where each genome sits in it, CanSNPCreator finds those
substitutions, reports which of them sit in well behaved sequence, and exports
a scheme for
[CanSNPer2-data](https://github.com/FOI-Bioinformatics/CanSNPer2-data).

## Installation

Python 3.13 or later; the package has no dependencies.

```sh
pip install git+https://github.com/FOI-Bioinformatics/CanSNPCreator.git
```

The `trf` subcommand also needs Tandem Repeats Finder on the PATH
(`conda install -c bioconda trf`).

## File formats

Every file is in a standard format, one per kind of data, so that the files
work with other tools as well as with this one.

| Data | Format | Notes |
| ---- | ------ | ----- |
| reference genome | FASTA | One record per contig. A contig is named by its ID, the first word of its header, as GFF3, VCF and BED name it. |
| alignment | FASTA | The reference genome as the first record, every record the same length. Gapped, where the reference row carries the aligner's gaps, or projected onto reference coordinates, as x2fa, snippy and RepGenR2 write it. |
| hierarchy | `child`/`parent` table, or Newick | The scheme's tree. The table is tab separated with a header naming `child` and `parent`, as FlexTaxD's taxonomy and RepGenR2's `tree2tax.tsv`; other columns are ignored. In Newick every label is a node name. |
| genome map | two tab separated columns, no header | Genome and the node it sits at, as FlexTaxD's `--genomeid2taxid` and RepGenR2's `genomes_map.tsv`. A genome may sit at an internal node. |
| genes | GFF3, from Prokka, Bakta or NCBI | `gene` and `CDS` features say whether a SNP is in a gene; each CDS's phase gives the reading frame. Described below. |
| repeats, masked regions | BED, or GFF3 | BED is 0-based and half open. Regions are named by contig ID, or, with `--masks-concatenated`, given over the whole concatenated reference, as Gubbins reports recombinant regions for an alignment in reference coordinates such as RepGenR2's `full_alignment.fasta`. |
| SNPs and scores | VCF 4.3 | Described below. |
| CanSNPer2-data scheme | that repository's three tab separated files | Written by `export`, read by `audit` and `place`. |

A root given as its own parent, as NCBI's `nodes.dmp` gives it, is read as the
root.

The hierarchy and the genome map are kept apart, as NCBI Taxonomy keeps its
nodes apart from its accession table: a node may hold no genome yet, and a
genome may sit at an internal node.

**Positions** are reported by contig, 1-based, as VCF and published schemes
give them. A column of the alignment is not a position in the genome, so the
alignment is mapped onto the reference, and the tools refuse to run unless the
two agree base for base.

**Names.** Contigs are matched by ID exactly, so a GFF3 or BED file must use
the reference FASTA's IDs; one that uses none of them is refused. Genome
names are matched more leniently, since a genome is often written as a file
name in one place and a record ID in another: `strain-7` matches `strain_7`,
and an assembly accession such as `GCF_000007765.2` matches the one record
whose name contains it, as RepGenR2 names records
`Family_Genus_species_Accession` but lists accessions alone. Node names are
never altered, so `2.10` stays `2.10`.

### Gene annotations

Prokka's, Bakta's and NCBI's GFF3 are read alike, compressed or not. A CDS
gives the reading frame, from its 5' end after skipping the phase GFF3 gives,
on either strand; a CDS of a pseudogene (`pseudo=true`, as NCBI marks it) or
one using a genetic code whose codons translate differently from the standard
code (`transl_table` other than 1 or 11) still places a SNP in a gene but gives
no frame. Prokka writes CDS features without gene features, which is enough.

The annotation's sequence names are matched to the reference's contigs in
three ways, in order:

1. by ID, when the annotation is of the reference FASTA it was made from, as
   Bakta's and Prokka's usually are and NCBI's are for the same accession;
2. by `--contig-aliases`, for the same assembly under other names: NCBI's
   `assembly_report.txt`, which pairs each GenBank accession with its RefSeq
   one, or a UCSC `chromAlias` table;
3. by the sequence itself, when the GFF3 carries it after `##FASTA`, as Prokka
   and Bakta write it, so that contigs Prokka renamed (`--compliant`) are
   matched by being identical.

An annotation none of whose sequences matches a contig is refused.

### SNPs and scores as VCF

`extract` writes every SNP of the alignment and `evaluate` the scored canSNPs,
both as VCF 4.3:

- `REF` is the reference base and `ALT` the other allele.
- `INFO/AA` is the ancestral allele, a key the VCF specification reserves for
  it; the derived allele is the other one. This holds whether or not the
  reference genome itself carries the derived allele.
- `INFO/NODE` is the node whose genomes carry the derived allele. A SNP that
  marks no node has neither key; `extract` keeps it, because it still counts as
  a neighbour when the others are scored.
- `INFO/NEAREST` is the distance, in bases, to the nearest other column where
  the genomes placed in the scheme differ, capped at 1000. It counts every such
  column, also those the strict core leaves out (a third base, or a genome
  without a base), and leaves the outgroup out.
- `INFO/MISSING` is the number of genomes without a base at the SNP, which is
  0 unless `extract --min-called` is below 1.
- `evaluate` adds `INFO/SCORE` and the five terms (`NEIGHBOUR`, `GC`, `TANDEM`,
  `GENE`, `MUTATION`), and records the weights in a `##CanSNPCreatorWeights`
  line.
- `INFO/COPIES`, from `evaluate`, is the highest number of times a k-mer
  holding the SNP occurs in the reference, both strands counted: 1 in unique
  sequence, more in a repeat, 0 when every such k-mer holds a base other than
  A, C, G or T.
- A canSNP that fails is kept, with a `FILTER` code: `LowScore` (not above
  `--cutoff`), `NotTop` (not among the `--top-n` best of its node that lie
  `--min-spacing` apart) or `Masked` (in a region given with `--masks`).
- `##contig` lines give each contig's ID and length. A VCF read against
  another reference is refused, and so is one whose `REF` disagrees with the
  reference anywhere.

```sh
bcftools view -f PASS -i 'INFO/NODE="B.Br.004"' scores.vcf
```

## Making a new scheme

The hierarchy and the genome map below are written by `collapse` from a genome
tree, as [Starting from a genome tree](#starting-from-a-genome-tree) describes,
or by hand. Every subcommand overwrites its output files.

```sh
# The canSNPs of every node, and every other SNP
CanSNPCreator extract \
    --alignment alignment.fasta --reference reference.fasta \
    --hierarchy hierarchy.tsv --genome-map genome_map.tsv \
    -o snps.vcf

# Tandem repeats of the reference, as BED
CanSNPCreator trf --reference reference.fasta -o repeats.bed

# Score them
CanSNPCreator evaluate \
    --snps snps.vcf --reference reference.fasta \
    --genes genes.gff3 --repeats repeats.bed --masks recombination.bed \
    --hierarchy hierarchy.tsv -o scores.vcf

# Write the scheme
CanSNPCreator export \
    --scores scores.vcf --hierarchy hierarchy.tsv --reference reference.fasta \
    --outdir scheme --organism francisella_tularensis \
    --genome SCHUS4.2 --strain francisella --citation Svensson2009 \
    --genbank-id GCA_000008985.1 --refseq-id GCF_000008985.1 --assembly-name ASM898v1 \
    --genome-map genome_map.tsv --alignment alignment.fasta
```

`extract` uses only sites where every genome has one of two bases (a strict
core): a genome with no base at a site, as an accessory region or a gap gives,
leaves the site out, so a large set of draft genomes leaves fewer sites than a
small set of complete ones.

`--min-called FRACTION` relaxes this. A site is kept when at least that
fraction of the genomes have a base and those genomes carry exactly two bases.
A node is marked when at least the same fraction of its own genomes have a
base, all of them the derived one, and no genome with a base outside the node
carries it. `INFO/MISSING` gives the number of genomes without a base, and
`evaluate` and `export` take, within a node, SNPs with fewer missing genomes
before better-scoring ones with more. With `--min-called 0.9`, a node of two
genomes still needs both, so the relaxation helps mainly the larger nodes and
sites where a genome outside the node lacks a base. The default, 1.0, keeps
the strict core.

### Scores

Five terms, each 0 to 100 before weighting:

| Term | Meaning |
| ---- | ------- |
| `neighbour_score` | distance to the nearest variation among the scheme's genomes (`INFO/NEAREST`), or to a contig end, over `--flank` bases (default 50); variation under a primer or probe interferes with an assay |
| `gc_score` | how close the GC content of the 100 bases around the SNP is to `--gc-goal`; full marks within `--gc-tolerance` (default 0.1) of it |
| `tandem_score` | distance to the closest tandem repeat |
| `gene_score` | whether the SNP lies in a gene or CDS, where its context is better defined |
| `mutation_score` | in the frame of the CDS it lies in: 100 if synonymous, 50 if it changes the amino acid, 0 if a stop codon is gained or lost or if it lies outside every CDS |

The gene term counts twice by default, so the highest total is 600;
`--score-weight gene_score=1` changes that. `--cutoff` defaults to 500, which in
practice requires the SNP to lie in a gene. The run reports the nodes left
without a passing canSNP, and those with fewer than the three CanSNPer2 needs
on a path to call a node.

A SNP file without `INFO/NEAREST`, from an earlier version, is scored by the
distance to the other SNPs of the file, which only sees the strict core. In a
set of 91 *Brucella* genomes that missed variation within 25 bases of more
than half of the canSNPs.

`evaluate` also counts how often the sequence around each canSNP occurs in the
reference: every k-mer of `--copy-kmer` bases (default 31) that holds the SNP
is looked up on both strands, and `INFO/COPIES` gives the highest count. A
canSNP in repeated sequence, such as an rRNA operon, an insertion sequence or a
paralogous gene, is harder to align to the right copy in a query genome. In
1,229 *Brucella* genomes left out of a scheme's tree, the 0.2% of candidates in
repeated 31-mers lacked a base three times as often as the others, and
disagreed with the other canSNPs of their node three times as often. Such
canSNPs are kept, and ranked after those in unique sequence. `--copy-kmer 0`
leaves the count out. Counting 187,000 candidates on a 3.4 Mbp reference takes
about five seconds and 0.7 GB of memory.

Within a node, canSNPs are taken in this order: unique sequence first, then
fewer genomes without a base (`INFO/MISSING`), then the higher score, then by
contig and position. In the same *Brucella* genomes, ranking further by how
consistently each candidate behaved in the held-out genomes changed the rank-1
SNP of few nodes and did not lower how often it disagreed with its node, so
this order is kept.

`--top-n` passes the best N canSNPs of each node, in that order, passing over
any within `--min-spacing` bases (default 1000) of one already taken, so that
the canSNPs of a node do not share a locus that one deletion or recombination
event could remove.

### Export

`export` takes the highest scoring passing canSNP of each node and writes
`<organism>_tree.txt`, `<organism>_snp.txt` and `<organism>_references.txt` as
CanSNPer2-data expects them, with positions counted over the concatenated
genome, since its SNP file has no contig column. That repository's rules are
checked first, so the result passes its `validate.py` unchanged. A node with no
canSNP refuses the scheme; `--prune` leaves such nodes out, joining their
children to the nearest node that has one.

CanSNPer2 calls one SNP per node, so the SNP file always holds the best of
each. `--snps-per-node N` (or `all`) keeps up to N canSNPs of each node, chosen
as `evaluate --top-n` chooses them, at least `--min-spacing` bases apart, and
writes them all to `<organism>_canSNPs.txt`:

```
node	rank	genome	position	derived_base	ancestral_base	score
```

ranked within each node, with rank 1 the SNP of the SNP file. Positions count
over the concatenated genome, as in the SNP file. CanSNPer2-data's
`validate.py` checks that rank 1 agrees with the SNP file, and `audit` checks
every rank against the reference and the genomes.

With `--genome-map` and `--alignment` it also writes `<organism>_genomes.txt`,
the genomes the scheme was made from: each by its versioned assembly
accession, the node it sits at, and its role, `representative` (its sequence
was in the alignment), `member` (placed through a representative, as
dereplication folds it in) or `outgroup` (aligned, at no node). A genome whose
name holds no accession is left out with a warning. CanSNPer2 does not read
the file; CanSNPer2-data's `validate.py` checks it, and `audit` uses it, so
that a published scheme carries what is needed to check it again.

## Starting from a genome tree

A genome tree, such as RepGenR2's `tree/tree.nwk`, has a node for every split
the tree builder found: far more than a scheme has, and many of them weakly
supported. `collapse` removes the internal nodes below a threshold, joining
each into its parent, and places every genome at the deepest node that
remains:

```sh
CanSNPCreator collapse \
    --tree tree.nwk --outgroup Coxiella_sp_GCF_000000000.1 \
    --min-support 0.95 --min-genomes 3 \
    --members genomes_map.tsv --rename Cb --outdir draft --prefix coxiella
```

| Option | Removes an internal node when |
| ------ | ----------------------------- |
| `--min-support` | its support is lower, as a fraction; a tree written in percentages is read as fractions |
| `--min-branch-length` | the branch above it is shorter, in the tree's units |
| `--min-genomes` | fewer genomes lie below it (default 2) |

A node holding exactly its parent's genomes is joined into the parent too.
`--outgroup` names a leaf placed at no node, so that it stays in the alignment
as a genome outside every node: that is what lets `extract` find the root's
canSNP and tell which side of each split is derived. Without one, `extract`
warns that the two sides of the root cannot be told apart. `--members` reads
RepGenR2's `genomes_map.tsv`, so that genomes dereplication folded into a
representative count towards its node and are placed with it. `--rename Cb`
names the nodes `Cb`, `Cb.1`, `Cb.1.2` and so on, largest first.

A Newick genome tree carries its support values and branch lengths; its
numeric internal labels are read as support. For `tree2tax.tsv`,
`--nodes-file tree2tax_nodes.tsv` supplies them. RepGenR2 keeps the outgroup as
a leaf in both unless `tree2tax --remove-outgroup` was given, and names it in
`outgroup_accession.txt`.

It writes `PREFIX_hierarchy.tsv` and `PREFIX_genomes.tsv`, the `--hierarchy`
and `--genome-map` of `extract`, and `PREFIX_nodes.tsv`, with each node's
genome count, support, branch length and name in the tree.

## Checking an existing scheme

### Audit

```sh
CanSNPCreator audit \
    --scheme-dir CanSNPer2-data/source --organism francisella_tularensis \
    --reference SCHUS4.2.fasta \
    --alignment alignment.fasta --tree tree.nwk \
    -o findings.tsv
```

CanSNPer2-data's `validate.py` answers whether a scheme is well formed; the
audit answers whether it can be typed with, and how it would be extended or
modified.

| Check | Needs | What it reports |
| ----- | ----- | --------------- |
| tree | | a node under two parents, a node with no SNP or two, a SNP naming a node the tree does not hold |
| alleles | | a base outside A, C, G and T, which can never be called derived; equal derived and ancestral bases; two nodes at one position |
| references | | a genome name with an underscore, a malformed accession, a SNP on a genome the references file lacks |
| callability | | nodes with fewer than `--min-required-hits` SNPs on the path from the root, which CanSNPer2 cannot call |
| spacing | | two canSNPs closer than `--proximity` bases, where one falls inside the other's assay |
| reference base | `--reference` | a position whose reference base is neither allele; when most are, coordinates taken against another version of the assembly |
| masked | `--reference`, `--masks` | a SNP in a masked region, such as a recombinant one |
| better SNPs | `--reference`, `--scores` | a node where a higher scoring passing canSNP is available, or one the scores do not cover |
| genomes | `--reference`, `--alignment`, and `--genome-map`, `--tree` or the scheme's own genomes file | the checks below |

When the scheme has an `<organism>_canSNPs.txt`, every rank in it is checked
by the reference base, masked, spacing and genome checks, and a finding about
a rank other than 1 says which rank it is about. Only a failing rank 1, the SNP
CanSNPer2 calls, is offered a replacement.

The genome checks are the only ones that test what defines a canSNP. For each
SNP of the scheme the audit reports genomes below the node without the derived
base (`not-canonical`) and genomes outside it with it (`homoplasy`), by name
and by the node they sit at; swapped alleles (`polarity`); genomes with no base
at the site (`missing-data`); for a failed SNP, other sites that separate
exactly that node's genomes (`replacement`), or that none does
(`unsupported-node`); a node that one site alone supports (`single-site`); and
a SNP that every genome carries, which nothing tests without an outgroup
(`no-outgroup`).

Each genome is also typed by its own canSNPs. A node is derived in a genome
when more of its canSNPs carry the derived base than the ancestral one, and the
genome's path ends at the deepest derived node with no ancestral node above it.
Two findings come from this. `misplaced` lists genomes whose canSNPs agree with
each other but lead to another node than the genome map gives, grouped by both
nodes. A genome map that places every member of a RepGenR2 cluster at its
representative's node gives many of these when the clusters are coarser than
the scheme: in a *Brucella* scheme of 91 nodes, 941 of 1,229 cluster members
were typed at another node than their representative: 345 above it, 234
below it and 362 on another branch. The finding gives the number of canSNPs the typing rests on, since
with one canSNP at a node the same pattern comes from a homoplastic SNP.
`discordant-genome` lists genomes that contradict the scheme at
`--genome-conflicts` nodes or more (default 3): derived off their path, or
ancestral on it. Base errors, contamination or a mixed assembly show
themselves this way, and such genomes are weak evidence against a SNP. The
`not-canonical` and `homoplasy` findings still count every genome, so a
`misplaced` or `discordant-genome` finding can explain them.

It also reads every variable site at which all genomes have one of two bases,
and reports groups of genomes that share SNPs no other genome carries and that
the scheme has no node for (`new-clade`, with the node they would hang from and
how many SNPs mark them), and sites that fit no grouping the scheme can hold
(`conflicting-sites`). A grouping that contradicts a better supported one is
counted as conflicting rather than offered. On a 2 Mbp genome, eight genomes
take about five seconds; the time grows with the number of genomes.

Each finding is a problem, which makes part of the scheme wrong or impossible
to type, or a recommendation, which leaves it usable, and each asks to
`modify` an existing SNP, `extend` the scheme, or `check` something the tool
cannot settle. A problem exits 2; a recommendation exits 0 unless `--strict` is
given, for continuous integration. `-o` writes the findings as a tab separated
file as well.

### Placing genomes by the tree

The genome checks need each genome placed at a node of the scheme. Typing the
genomes with the scheme places them by the SNPs under test, so a genome lacking
the derived base of its own branch is typed elsewhere and never seen to lack
it. `--tree`, like the `place` subcommand, uses each SNP only to identify which
clade of the genome tree the node means, and places the genomes by the tree:

```sh
CanSNPCreator place \
    --scheme-dir CanSNPer2-data/source --organism coxiella_burnetii \
    --tree tree.nwk --alignment alignment.fasta --reference RSA493.fasta \
    -o placed.tsv --report placement.tsv
```

Each node is matched to the side of a branch whose genomes best agree with the
genomes carrying its derived base, by Jaccard index, counting only genomes
with a base at the site. Both sides of every branch are considered, so where
the tree is rooted does not matter. A genome is placed at the deepest node
whose match holds it. `--min-match` (default 0.5) sets the lowest index
accepted, and the report gives each node's match and index. A node with one or
two genomes is matched less reliably than one with many: in a clade of two in
which one genome lacks the derived base, the other alone matches as well as the
clade.

## Exit status

| Status | Meaning |
| ------ | ------- |
| 0 | finished, output written |
| 1 | a file could not be read or written, an input is not usable, or options were combined that do not go together |
| 2 | the inputs were read but are inconsistent with one another, the audited scheme holds a problem, `place` placed no genome, or no subcommand was given; also an option argparse refuses, such as a value out of range |

`--log-level DEBUG` reports more, `--log-level WARNING` only what went wrong;
it is given before the subcommand.

## Development

```sh
pip install -e ".[dev]"
python -m pytest tests
python -m ruff check cansnpscorer tests
```

The tests build a small dataset with one canSNP planted per node, including a
two-contig reference, a node whose genomes are two, a SNP crossing two
branches, a tri-allelic column and a gapped one, and check that the planted
canSNPs are recovered at the planted coordinates with the planted polarity. CI
runs them on Python 3.13 and 3.14.

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
