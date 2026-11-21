# Calling a node from several canSNPs

This document describes how CanSNPCreator types a genome when a scheme holds
more than one canSNP per node (the "set rule"), and what a CanSNPer2 call from a
set would involve. The rules below are those implemented in
`cansnpscorer/typing_report.py`; the function that implements each is named. It
describes existing behaviour and makes no change to CanSNPer2.

## 1. Input

`export --snps-per-node N` writes `<organism>_canSNPs.txt` next to the usual
`<organism>_tree.txt` and `<organism>_snp.txt` (README, section "Export"):

```
node	rank	genome	sequence	sequence_position	derived_base	ancestral_base	score
```

Rank is counted within each node. Rank 1 is the SNP of `<organism>_snp.txt`, so a
reader that ignores the set sees the single-SNP scheme. Each canSNP is given by
its reference sequence and a 1-based position on it, as in the SNP file
(CanSNPer2 2.3.0). With the default
`--snps-per-node 1` the file holds rank 1 only.

## 2. State of a node (`snp_state`, `node_state`)

Each canSNP of a genome has one of four states from the base at its column:
derived, ancestral, other (a base that is neither allele) or missing (no base,
gap or ambiguity code).

The state of a node is the majority of its canSNPs that are derived or ancestral:

- more derived than ancestral: derived;
- more ancestral than derived: ancestral;
- equal, or none of either: missing.

A canSNP in the state "other" or "missing" does not vote. Three canSNPs
read as derived, missing, missing give a derived node; one derived and one
ancestral give a missing node.

## 3. The call (`set_path`, `set_call`)

1. A node is a candidate when it is derived and no node above it (towards the
   root) is ancestral. A node of missing state does not block a candidate.
2. For each candidate, the support is the number of derived nodes on its lineage,
   including itself.
3. The candidates with the greatest support are kept. If one of them is deeper
   than, and below, all the others (all others lie on its lineage), it is the
   call. Otherwise the best candidates lie on different branches and there is no
   call (reported as NA). The path is then the part of their lineages they share.
4. With no candidate there is no call and the path is empty.

The set rule has no minimum chain length: a single derived root node can be a
call. This differs from the chain requirement of the single rule below.

The call is not the deepest node on a chain of derived nodes in the strict
sense. A missing node on the way is passed, and the choice is by number of
derived nodes first and depth second.

For comparison, the single rule (`single_call`) reads only the rank 1 canSNP of
each node, counts a node as derived only when its rank 1 canSNP is derived, and
follows CanSNPer2 2.2.1 and 2.2.2 (section 4, `confirm_chain`): an ancestral or third-base
rank 1 canSNP on the lineage breaks the chain, a missing one is a gap, a node
without a SNP is passed and not counted, and the
call is the deepest derived node with at least `min_required_hits` derived nodes
on its lineage (the root counts as one), a derived fraction of at least
`strictness`, and a derived parent. Two deepest confirmed nodes at the same depth
give no call. The single-rule figures in section 6 were measured with the rule of
CanSNPer2 2.2.0, where a missing rank 1 SNP ended the chain.

## 4. `--min_required_hits` and `--strictness` under the set rule

CanSNPer2 calls one SNP per node. The rule below is that of CanSNPer2 2.2.1
(commit 2e8391c, PR #55), `NewickTree._confirm_path`, and is unchanged in
2.2.2 (upstream master, commit dc660a9, PRs #59 and #60). Each SNP
on the chain from the root to the node is one of:

- called: derived in the query, counted as a hit;
- no data: no alignment covered the position, or the query has a gap or N there.
  In 2.2.2 a position that two alignments cover with different query bases,
  which can occur with `--aligner minimap2`, is also read as N. This is a gap in
  the path, not a break;
- against: ancestral or a third base. The chain is broken and the node cannot be
  called.

A node is confirmed when the chain is not broken, at least `--min_required_hits`
(default 3) SNPs are called, the called fraction of the SNPs on the chain is at
least `--strictness` (default 0.7), and the SNP directly above the node is
called.

Under the set rule, the corresponding reading in CanSNPCreator is:

- The unit of the minimum is a node, as before, but a node is derived from the
  majority of its canSNPs. `min_required_hits` and `strictness` are applied by
  the single rule (`single_call`, `audit --min-required-hits --strictness`) and
  not by the set rule.
- A node without a majority (missing state) is passed by the set rule. This is
  the same kind of tolerance as `--strictness` gives for nodes without alignment
  data, and a node whose canSNPs are all missing is treated alike. A node in the
  ancestral state is not passed.
- The set rule counts the missing nodes on the path (`missing_on_path`), so that
  a call resting on many unknown nodes can be recognised.

In 2.2.1 and 2.2.2, `--min_required_hits` is a minimum count of called SNPs and
`--strictness` is the tolerated fraction of SNPs without data on the chain; this
already covers a path with unknown nodes. What a set rule would add is the state
of a node from several SNPs (section 2), so that one missing or discordant base
no longer makes the node a gap or a break, and the per-node counts of section 5.
CanSNPCreator's set rule has no counterpart of the count and fraction thresholds
and does not require the node above the call to be derived.

## 5. Output per node that CanSNPer2 would add

From the set rule, `type_alignment` records for each genome the following, and
`audit` writes them (`<prefix>_genomes.tsv`):

- `derived_on_path`: nodes on the path that are derived;
- `ancestral_on_path`: nodes on the path that are ancestral (cannot occur on a
  path that has a call, since an ancestral node blocks it; kept as a check);
- `missing_on_path`: nodes on the path that are missing, including nodes with no
  canSNP in the scheme;
- `derived_off_path` (`derived_off`): derived nodes that contradict the path. With
  a call these are the derived nodes off its path. Where branches tie and there is
  no call, it is the fewest that contradict any one of the tied branches.

Per node, a CanSNPer2 report would need the number of its canSNPs that were
derived, ancestral and missing in the genome, so that the majority and its margin
can be read. CanSNPCreator reports the same counts per canSNP and node in
`<prefix>_snps.tsv` and `<prefix>_nodes.tsv`.

## 6. Benchmark evidence

Held-out genomes were projected onto each organism's alignment and typed by both
rules with `--snps-per-node 5` (CanSNPCreator `audit`), the single rule being
that of CanSNPer2 2.2.0 at the time. The table gives the
number of held-out genomes, the number typed (a call) by each rule where stated,
and the share of typed genomes where the single call equals the set call. The
table lists the organisms with a single = set agreement figure in the benchmark
notes, plus Brucella for its typed counts. Figures are as recorded in the
benchmark notes. The Brucella typed counts are over an alignment of the
reference plus the 1,229 held-out genomes, which is why the set rule count
(1,230) exceeds the held-out count. "-" means not stated.

| Organism | Held out | Single rule typed | Set rule typed | Single = set |
|---|---|---|---|---|
| Legionella pneumophila | 834 | 833 | - | 95.1% |
| Listeria monocytogenes | 920 | 916 | - | 95.5% |
| Mycobacterium tuberculosis | 1,113 | 1,112 | - | 95.2% |
| Vibrio cholerae | 819 | 817 | - | 98.8% |
| Staphylococcus aureus | 767 | 766 | - | 94.4% |
| Helicobacter pylori | 174 | 174 | 174 | 99.4% |
| Yersinia pestis (159 representatives) | 880 | 873 | 880 | 91.1% |
| Yersinia pestis (386 representatives) | 653 | 651 | 653 | 91.6% |
| Coxiella burnetii | 127 | 124 | - | 92.1% |
| Escherichia coli | 279 | 279 | - | 92.1% |
| Shigella sonnei | 1,314 | 1,314 | - | 92.9% |
| Salmonella Typhi | 1,033 | - | - | 84.8% |
| Brucella melitensis | 1,229 | 1,224 | 1,230 | - |

Where the two calls differ, the set call is usually deeper: for Brucella 39
genomes were typed deeper by the set rule (21 where the rank 1 SNP was ancestral,
15 missing, 3 for `min_required_hits`). In Y. pestis the single rule stopped at a
canSNP without a base in 70 genomes (8.0%), and in S. Typhi at a missing base at
one node for 31 genomes. In the recorded notes, the K. pneumoniae genomes not
typed by the single rule (23 of 266) all stopped for too few nodes on the chain
(`min_required_hits` 3).

The set call is not an independent truth. It reads the same genomes and the same
alignment as the single call, and the table shows how often the two agree, not
which is correct. Agreement with the node expected from the dereplication
clusters was lower and depends on whether the clusters follow the clades. After
`evaluate --observed` (split-half test), agreement of single with set rose in the
test half for most organisms (for example 95.9% to 98.4% in Brucella), with
exceptions recorded in the notes.
