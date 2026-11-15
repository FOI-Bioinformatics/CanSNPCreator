import argparse
import logging
import sys
import zlib

from cansnpscorer import (
    __version__,
    cansnp_audit,
    cansnp_evaluator,
    cansnp_export,
    cansnp_extractor,
    cansnp_tree,
    cansnp_trf,
)
from cansnpscorer.alignment_audit import GENOME_CONFLICTS
from cansnpscorer.helpers.data_utils import MINIMUM_SNPS_PER_NODE
from cansnpscorer.helpers.kmers import DEFAULT_K
from cansnpscorer.helpers.scoring_methods import (
    DEFAULT_PROFILE,
    MAXIMUM_SCORE,
    PROFILES,
    TERMS,
)


def at_least(minimum: int):
    """An argparse type for a whole number no smaller than minimum."""
    def parse(text: str) -> int:
        try:
            value = int(text)
        except ValueError:
            raise argparse.ArgumentTypeError(f"'{text}' is not a whole number") from None
        if value < minimum:
            raise argparse.ArgumentTypeError(f"must be at least {minimum}, not {value}")
        return value
    return parse


def count_or_all(text: str):
    """An argparse type for a positive whole number, or 'all' for None."""
    if text.lower() == "all":
        return None
    return at_least(1)(text)


def between(low: int, high: int):
    """An argparse type for a whole number from low to high."""
    def parse(text: str) -> int:
        value = at_least(low)(text)
        if value > high:
            raise argparse.ArgumentTypeError(f"must be at most {high}, not {value}")
        return value
    return parse


def fraction(text: str) -> float:
    """An argparse type for a number from 0 to 1."""
    try:
        value = float(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"'{text}' is not a number") from None
    if not 0.0 <= value <= 1.0:
        raise argparse.ArgumentTypeError(f"must be from 0 to 1, not {value:g}")
    return value


def called_fraction(text: str) -> float:
    """A fraction above 0 and at most 1."""
    value = fraction(text)
    if value <= 0:
        raise argparse.ArgumentTypeError("must be above 0")
    return value


def non_negative_float(text: str) -> float:
    try:
        value = float(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"'{text}' is not a number") from None
    if value < 0:
        raise argparse.ArgumentTypeError(f"must not be negative, not {value:g}")
    return value


def score_weights(settings, profile=DEFAULT_PROFILE):
    """Turn --score-weight TERM=FACTOR arguments into a dictionary of weights.

    The weights are those of the profile, with the terms the settings name
    changed; a term not named keeps the profile's weight.
    """
    weights = dict(PROFILES[profile][0])
    if not settings:
        return weights

    for setting in settings:
        term, _, factor = setting.partition("=")
        if term not in TERMS:
            raise argparse.ArgumentTypeError(
                f"unknown score term '{term}'; choose from {', '.join(TERMS)}")
        try:
            weights[term] = float(factor)
        except ValueError:
            raise argparse.ArgumentTypeError(
                f"the weight for '{term}' must be a number, got '{factor}'") from None
    return weights


def add_version(parser) -> None:
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")


def add_reference(parser, required: bool = True, purpose: str = "") -> None:
    parser.add_argument(
        "--reference", required=required, default=None, metavar="FASTA",
        help="The reference genome in FASTA format, one record per contig, in the order the "
             "alignment concatenates them. A contig is named by its ID, the first word of its "
             "header, as GFF3, VCF and BED name it." + (f" {purpose}" if purpose else ""))


def add_alignment(parser, required: bool = True) -> None:
    parser.add_argument(
        "--alignment", required=required, default=None, metavar="FASTA",
        help="The whole genome alignment in FASTA format, with the reference as its first "
             "record: either gapped, or projected onto reference coordinates as x2fa, snippy "
             "and RepGenR2 write it.")


def add_hierarchy(parser, required: bool = True, purpose: str = "") -> None:
    parser.add_argument(
        "--hierarchy", required=required, default=None, metavar="FILE",
        help="The scheme's tree: a tab separated table with 'child' and 'parent' columns, as "
             "FlexTaxD and RepGenR2's tree2tax.tsv, or Newick whose labels name the nodes."
             + (f" {purpose}" if purpose else ""))


def add_genome_map(parser, required: bool = True, purpose: str = "") -> None:
    parser.add_argument(
        "--genome-map", required=required, default=None, metavar="TSV",
        help="Which genome sits at which node: two tab separated columns, genome and node, "
             "without a header, as FlexTaxD's --genomeid2taxid. A genome may sit at an "
             "internal node, and may be named by its alignment record or its accession."
             + (f" {purpose}" if purpose else ""))


def add_members(parser) -> None:
    parser.add_argument(
        "--members", default=None, metavar="TSV",
        help="RepGenR2's genomes_map.tsv, which maps each genome dereplication folded into a "
             "representative to that representative's leaf, so that it is placed with it.")


def add_masks(parser, purpose: str) -> None:
    parser.add_argument(
        "--masks", default=None, metavar="BED",
        help="Regions to mask, such as recombinant ones, in BED or GFF3 named by contig ID. "
             + purpose)
    parser.add_argument(
        "--masks-concatenated", action="store_true",
        help="Read the positions of --masks over the whole concatenated reference, whatever "
             "sequence name the file gives, as Gubbins reports them for an alignment in "
             "reference coordinates such as RepGenR2's full_alignment.fasta.")


def add_contig_aliases(parser) -> None:
    parser.add_argument(
        "--contig-aliases", default=None, metavar="FILE",
        help="Other names of the reference's contigs, so that annotation files from another "
             "naming of the same assembly can be used: NCBI's assembly_report.txt, which pairs "
             "GenBank and RefSeq accessions, or a UCSC chromAlias table.")


def main():
    parser = argparse.ArgumentParser(
        prog="CanSNPCreator",
        description="Find, score, export and audit canonical SNPs from whole genome alignments.")
    add_version(parser)
    parser.add_argument("--log-level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        help="How much to report while running; given before the subcommand. "
                             "Default: INFO.")
    subparsers = parser.add_subparsers(dest="command")

    # extract
    extract = subparsers.add_parser(
        "extract", help="Find the SNPs of an alignment, and the node each canonical one marks.")
    add_alignment(extract)
    add_reference(extract)
    add_hierarchy(extract)
    add_genome_map(extract)
    extract.add_argument("--min-called", type=called_fraction, default=1.0, metavar="FRACTION",
                         help="The fraction of genomes, overall and among a node's genomes, "
                              "that must have a base (not a gap or N) at a SNP. Default 1.0: "
                              "a column with any missing base is left out. Below 1.0 the "
                              "genomes without a base are ignored when the node is judged, "
                              "INFO/MISSING counts them, and export prefers SNPs with fewer.")
    extract.add_argument("-o", "--output", required=True, help="The SNPs to write.")
    add_version(extract)

    # evaluate
    evaluate = subparsers.add_parser(
        "evaluate", help="Score and rank canSNPs for how reliably they can be called.")
    evaluate.add_argument("--snps", required=True, help="The SNPs extract wrote.")
    add_reference(evaluate)
    evaluate.add_argument("--genes", required=True, metavar="GFF3",
                          help="The reference's genes in GFF3, as Prokka, Bakta or NCBI write "
                               "them; CDS features give the reading frame, gene and CDS features "
                               "say whether a SNP is in a gene. Contigs Prokka or Bakta renamed "
                               "are matched by the sequences they carry after ##FASTA.")
    evaluate.add_argument("--repeats", required=True, metavar="BED",
                          help="The reference's tandem repeats, in BED as the trf subcommand "
                               "writes them, or in GFF3.")
    add_masks(evaluate, "A canSNP in one is kept with FILTER=Masked.")
    add_contig_aliases(evaluate)
    add_hierarchy(evaluate, required=False,
                  purpose="With it, the nodes left without a canSNP are reported.")
    evaluate.add_argument("-o", "--output", required=True, help="The scored canSNPs to write.")
    evaluate.add_argument("--tandem-limit", type=at_least(0), default=50,
                          help="How near a tandem repeat still lowers the score, in bases. "
                               "Default: 50.")
    evaluate.add_argument("--gene-limit", type=at_least(0), default=0,
                          help="How near a gene still raises the score, in bases. Default: 0.")
    evaluate.add_argument("--gc-goal", type=fraction, default=0.5,
                          help="The GC content that scores highest. Default: 0.5.")
    evaluate.add_argument("--gc-tolerance", type=fraction, default=0.1,
                          help="A GC content within this of --gc-goal scores as the goal "
                               "itself. Default: 0.1.")
    evaluate.add_argument("--flank", type=at_least(1), default=50,
                          help="How many bases on either side of a canSNP must be free of "
                               "variation among the scheme's genomes, and of a contig end, for "
                               "the full neighbour term. Default: 50. Distances beyond 1000 "
                               "bases are not recorded by extract.")
    evaluate.add_argument("--profile", choices=sorted(PROFILES), default=DEFAULT_PROFILE,
                          help="What the scheme is for, which sets the weights and the cutoff "
                               "used unless --score-weight and --cutoff give others. 'typing', "
                               "for schemes typed from genome alignments: no GC term, the gene "
                               "term once, and no cutoff, so that every node keeps the canSNPs "
                               "it has and they are only ranked. 'assay', for schemes typed by "
                               "PCR assays, and the behaviour of versions up to 0.3: every "
                               f"term, the gene term twice, cutoff 500. Default: "
                               f"{DEFAULT_PROFILE}.")
    evaluate.add_argument("--cutoff", type=at_least(0), default=None,
                          help=f"Only pass canSNPs scoring above this. The highest reachable "
                               f"total is {MAXIMUM_SCORE} with the assay profile's weights, "
                               f"whose default of 500 in practice requires the SNP to fall in "
                               f"a gene. Default: that of --profile.")
    evaluate.add_argument("--observed", default=None, metavar="TSV",
                          help="The _snps.tsv table audit --typing-report wrote for a scheme "
                               "exported from these canSNPs. CanSNPs that disagreed with their "
                               "node, or had no base, in the genomes typed are ranked after "
                               "the others of their node, ahead of the score. Where the "
                               "canSNPs of a node are derived in groups of genomes of "
                               "different width and no group holds more than half of them, "
                               "they are ranked by their distance to the narrowest group.")
    evaluate.add_argument("--top-n", type=at_least(1), default=None,
                          help="Only pass the best this many canSNPs of each node.")
    evaluate.add_argument("--min-spacing", type=at_least(0), default=1000,
                          help="With --top-n, the fewest bases between two canSNPs passed for "
                               "the same node. Beyond it, canSNPs in different genes are taken "
                               "before ones sharing a gene. Default: 1000.")
    evaluate.add_argument("--copy-kmer", type=at_least(0), default=DEFAULT_K, metavar="K",
                          help=f"Count how often the sequence around each canSNP occurs in the "
                               f"reference, with k-mers of this size on both strands, and rank "
                               f"canSNPs in repeated sequence after the others of their node. "
                               f"0 turns the count off. Default: {DEFAULT_K}.")
    evaluate.add_argument("--score-weight", action="append", metavar="TERM=FACTOR",
                          help=f"Weight one score term before the terms are summed. Terms: "
                               f"{', '.join(TERMS)}. Repeatable. Terms not named keep the "
                               f"weight of --profile: "
                               + "; ".join(
                                   f"{name}: " + ", ".join(
                                       f"{term}={weights.get(term, 1.0):g}" for term in TERMS)
                                   for name, (weights, _cutoff) in sorted(PROFILES.items()))
                               + ".")
    add_version(evaluate)

    # export
    export = subparsers.add_parser(
        "export", help="Write scored canSNPs as the source files of a CanSNPer2-data scheme.")
    export.add_argument("--scores", required=True, help="The scored canSNPs evaluate wrote.")
    add_hierarchy(export)
    add_reference(export, purpose="Its contig lengths turn positions into the concatenated "
                                  "coordinates the scheme's SNP file uses.")
    export.add_argument("--organism", required=True,
                        help="Names the three files, as in francisella_tularensis.")
    export.add_argument("--genome", required=True,
                        help="The short name of the reference genome, as the SNP file refers "
                             "to it. It may not hold an underscore.")
    export.add_argument("--strain", required=True, help="The strain column of both files.")
    export.add_argument("--citation", required=True,
                        help="What the SNPs are cited as, the reference column of the SNP file, "
                             "as in Svensson2009.")
    export.add_argument("--genbank-id", required=True, help="As in GCA_000007125.1.")
    export.add_argument("--refseq-id", required=True, help="As in GCF_000007125.1.")
    export.add_argument("--assembly-name", required=True, help="As in ASM712v1.")
    export.add_argument("--outdir", default=".", help="Where to write the files.")
    export.add_argument("--genome-map", default=None, metavar="TSV",
                        help="Where each genome sits, as extract used it. With --alignment, "
                             "the genomes the scheme was made from are written to "
                             "<organism>_genomes.txt, by assembly accession and role.")
    export.add_argument("--alignment", default=None, metavar="FASTA",
                        help="The alignment the canSNPs were found in, which tells the genomes "
                             "that were aligned from those placed through them.")
    export.add_argument("--snps-per-node", type=count_or_all, default=1, metavar="N|all",
                        help="How many canSNPs of each node to keep, best first. The SNP file "
                             "keeps the best, the one CanSNPer2 calls; with more than one, all of "
                             "them go to <organism>_canSNPs.txt, ranked within the node. "
                             "Default: 1.")
    export.add_argument("--min-spacing", type=at_least(0), default=1000,
                        help="The fewest bases between two canSNPs kept for the same node. "
                             "Beyond it, canSNPs in different genes are kept before ones "
                             "sharing a gene. Default: 1000.")
    export.add_argument("--prune", action="store_true",
                        help="Leave out nodes that have no canSNP, instead of refusing to write "
                             "a scheme. Their children then hang from the nearest node that has "
                             "one.")
    add_version(export)

    # audit
    audit = subparsers.add_parser(
        "audit", help="Report whether an existing CanSNPer2-data scheme should be extended or "
                      "modified.")
    audit.add_argument("--scheme-dir", required=True,
                       help="The directory holding the scheme's three source files, as the "
                            "source directory of CanSNPer2-data.")
    audit.add_argument("--organism", default=None,
                       help="Which scheme to read; may be left out when the directory holds "
                            "one.")
    add_reference(audit, required=False,
                  purpose="With it, the base at every position is checked against the two "
                          "alleles, the direct test for positions taken against another "
                          "version of the assembly.")
    audit.add_argument("--scores", default=None,
                       help="Scored canSNPs from evaluate, to report where a better one is "
                            "available. Needs --reference.")
    audit.add_argument("--genome", default=None,
                       help="Which of the scheme's reference genomes --reference holds, when "
                            "the scheme names more than one.")
    add_masks(audit, "SNPs of the scheme in them are reported. Needs --reference.")
    add_contig_aliases(audit)
    audit.add_argument("--min-required-hits", type=at_least(1), default=MINIMUM_SNPS_PER_NODE,
                       help=f"The SNPs CanSNPer2 requires on the path to a node before it calls "
                            f"it. Default: {MINIMUM_SNPS_PER_NODE}, as CanSNPer2 uses.")
    audit.add_argument("--proximity", type=at_least(0), default=50,
                       help="Report two canSNPs closer than this many bases. Default: 50.")
    audit.add_argument("--set-region", type=at_least(0), default=10000,
                       help="Report a node more than half of whose ranked canSNPs, and at least "
                            "three, lie within this many bases of each other; 0 turns the "
                            "check off. Default: 10000.")
    add_alignment(audit, required=False)
    add_genome_map(audit, required=False,
                   purpose="With --alignment, each SNP is checked in every genome. Without it "
                           "or --tree, the scheme's own genomes file is used.")
    audit.add_argument("--tree", default=None,
                       help="Instead of --genome-map: a genome tree whose leaves are the "
                            "alignment's genomes, which places them on the scheme by the tree, "
                            "as the place subcommand does, rather than by typing.")
    add_members(audit)
    audit.add_argument("--min-clade-size", type=at_least(1), default=2,
                       help="Report a clade the scheme does not resolve only when it holds at "
                            "least this many genomes. Default: 2.")
    audit.add_argument("--min-clade-snps", type=at_least(1), default=1,
                       help="Report a clade the scheme does not resolve only when at least this "
                            "many SNPs mark it. Default: 1.")
    audit.add_argument("--genome-conflicts", type=at_least(1), default=GENOME_CONFLICTS,
                       metavar="N",
                       help="With --alignment, report a genome whose canSNPs contradict the "
                            "scheme at this many nodes or more: derived off the path they type "
                            "it to, or ancestral on it. Default: "
                            f"{GENOME_CONFLICTS}.")
    audit.add_argument("--typing-report", default=None, metavar="PREFIX",
                       help="With --alignment, type every genome with the scheme, by the "
                            "rule of CanSNPer2 (one canSNP per node) and by the majority of "
                            "each node's canSNPs, and write PREFIX_genomes.tsv, _snps.tsv, "
                            "_nodes.tsv, _ranks.tsv, _summary.tsv and _run.tsv. No genome map "
                            "is needed; with one, "
                            "each call is compared with the node it gives. evaluate "
                            "--observed reads the _snps.tsv table.")
    audit.add_argument("--report", default=None, metavar="FILE",
                       help="Write the findings, and with --alignment the typing of its "
                            "genomes, as one Markdown report for reading.")
    audit.add_argument("-o", "--output", default=None,
                       help="Also write the findings to this tab separated file.")
    audit.add_argument("--strict", action="store_true",
                       help="Exit 2 for a recommendation as well as for a problem.")
    add_version(audit)

    # collapse
    collapse = subparsers.add_parser(
        "collapse", help="Reduce a genome tree to the nodes of a scheme, and place its genomes.")
    collapse.add_argument("--tree", required=True,
                          help="A genome tree whose leaves are genomes: Newick, or a child and "
                               "parent table such as RepGenR2's tree2tax.tsv.")
    collapse.add_argument("--nodes-file", default=None,
                          help="RepGenR2's tree2tax_nodes.tsv, for the support and branch length "
                               "of each node. A Newick tree carries both.")
    add_members(collapse)
    collapse.add_argument("--min-support", type=fraction, default=None,
                          help="Remove internal nodes with a lower support, as a fraction; a tree "
                               "in percentages is read as fractions.")
    collapse.add_argument("--min-genomes", type=at_least(1), default=2,
                          help="Remove internal nodes with fewer genomes below them. Default: 2.")
    collapse.add_argument("--min-branch-length", type=non_negative_float, default=None,
                          help="Remove internal nodes whose branch is shorter, in the tree's "
                               "units.")
    collapse.add_argument("--outgroup", action="append", default=[], metavar="LEAF",
                          help="A leaf that is an outgroup: it sits at no node and stays in the "
                               "alignment as a genome outside every node, which gives the root a "
                               "polarity. The leaf name, or an assembly accession it holds, as "
                               "RepGenR2's outgroup_accession.txt gives it. Repeatable.")
    collapse.add_argument("--rename", default=None, metavar="PREFIX",
                          help="Name the nodes PREFIX, PREFIX.1, PREFIX.1.2 and so on, largest "
                               "first and, among children of one size, by their first genome "
                               "name, instead of the tree's own names.")
    collapse.add_argument("--outdir", default=".", help="Where to write the files.")
    collapse.add_argument("--prefix", default="collapsed",
                          help="Names the files PREFIX_hierarchy.tsv, PREFIX_genomes.tsv and "
                               "PREFIX_nodes.tsv.")
    add_version(collapse)

    # place
    place = subparsers.add_parser(
        "place", help="Place the genomes of a tree on the nodes of an existing scheme.")
    place.add_argument("--scheme-dir", required=True,
                       help="The directory holding the scheme's three source files.")
    place.add_argument("--organism", default=None, help="Which scheme to read.")
    place.add_argument("--tree", required=True,
                       help="A genome tree whose leaves are the alignment's genomes; its root "
                            "does not matter.")
    add_alignment(place)
    add_reference(place)
    place.add_argument("--genome", default=None,
                       help="Which of the scheme's reference genomes --reference holds, when the "
                            "scheme names more than one.")
    add_members(place)
    place.add_argument("--min-match", type=fraction, default=0.5,
                       help="The lowest Jaccard index between a node's derived genomes and a "
                            "side of a branch that counts as a match. Default: 0.5.")
    place.add_argument("-o", "--output", required=True,
                       help="The genome map to write: genome and node, tab separated.")
    place.add_argument("--report", default=None,
                       help="Also write how each node was matched to the tree.")
    add_version(place)

    # trf
    trf = subparsers.add_parser("trf", help="Find the tandem repeats of the reference genome.")
    add_reference(trf)
    trf.add_argument("-o", "--output", required=True, metavar="BED",
                     help="The repeats to write, as BED.")
    trf.add_argument("--match", type=at_least(1), default=2, help="Matching weight. Default: 2.")
    trf.add_argument("--mismatch", type=at_least(1), default=7, help="Mismatch penalty. Default: 7.")
    trf.add_argument("--delta", type=at_least(1), default=7, help="Indel penalty. Default: 7.")
    trf.add_argument("--pm", type=between(0, 100), default=80, help="Match probability. Default: 80.")
    trf.add_argument("--pi", type=between(0, 100), default=10, help="Indel probability. Default: 10.")
    trf.add_argument("--min-score", type=at_least(1), default=10,
                     help="Minimum alignment score to report. Default: 10.")
    trf.add_argument("--max-period", type=between(1, 2000), default=2000,
                     help="Maximum period size to report. Default: 2000.")
    add_version(trf)

    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(asctime)s - %(levelname)s - %(message)s")
    if getattr(args, "masks_concatenated", False) and not getattr(args, "masks", None):
        parser.error("--masks-concatenated applies to --masks")
    try:
        return run(args, parser, evaluate)
    except (OSError, EOFError, zlib.error) as error:
        # A file that cannot be written, or a truncated or corrupt compressed
        # file, is reported as what it is rather than as a traceback.
        logging.error(f"{type(error).__name__}: {error}")
        return 1


def run(args, parser, evaluate) -> int:
    """Run the subcommand the arguments name. Returns an exit status."""

    if args.command == "extract":
        return cansnp_extractor.main(
            alignment=args.alignment, reference=args.reference, hierarchy=args.hierarchy,
            genome_map=args.genome_map, output=args.output, min_called=args.min_called)
    if args.command == "evaluate":
        try:
            weights = score_weights(args.score_weight, args.profile)
        except argparse.ArgumentTypeError as error:
            evaluate.error(str(error))
        cutoff = PROFILES[args.profile][1] if args.cutoff is None else args.cutoff
        reachable = round(100 * sum(weights.get(term, 1.0) for term in TERMS))
        if cutoff is not None and cutoff >= reachable:
            logging.warning(
                f"No canSNP can pass: a score must be above the cutoff of {cutoff}, and the "
                f"highest total these weights can give is {reachable}. Give a lower --cutoff")
        return cansnp_evaluator.main(
            snps=args.snps, reference=args.reference, genes=args.genes, repeats=args.repeats,
            hierarchy=args.hierarchy, output=args.output, tandem_limit=args.tandem_limit,
            gene_limit=args.gene_limit, gc_goal=args.gc_goal,
            cutoff=cutoff,
            top_n=args.top_n, weights=weights, masks=args.masks,
            masks_concatenated=args.masks_concatenated, contig_aliases=args.contig_aliases,
            flank=args.flank, gc_tolerance=args.gc_tolerance, spacing=args.min_spacing,
            copy_kmer=args.copy_kmer, observed=args.observed)
    if args.command == "export":
        return cansnp_export.main(
            scores=args.scores, hierarchy=args.hierarchy, reference=args.reference,
            outdir=args.outdir, organism=args.organism, genome=args.genome, strain=args.strain,
            citation=args.citation, genbank_id=args.genbank_id, refseq_id=args.refseq_id,
            assembly_name=args.assembly_name, prune=args.prune, genome_map=args.genome_map,
            alignment=args.alignment, snps_per_node=args.snps_per_node,
            spacing=args.min_spacing)
    if args.command == "audit":
        return cansnp_audit.main(
            scheme_dir=args.scheme_dir, organism=args.organism, reference=args.reference,
            scores=args.scores, genome=args.genome, min_required_hits=args.min_required_hits,
            proximity=args.proximity, output=args.output, strict=args.strict,
            alignment=args.alignment, genome_map=args.genome_map,
            min_clade_size=args.min_clade_size, min_clade_snps=args.min_clade_snps,
            tree=args.tree, members=args.members, masks=args.masks,
            masks_concatenated=args.masks_concatenated, contig_aliases=args.contig_aliases,
            genome_conflicts=args.genome_conflicts, typing_report=args.typing_report,
            report=args.report, region=args.set_region)
    if args.command == "collapse":
        return cansnp_tree.collapse_main(
            tree_file=args.tree, nodes_file=args.nodes_file, members=args.members,
            min_support=args.min_support, min_genomes=args.min_genomes,
            min_length=args.min_branch_length, prefix_names=args.rename, outdir=args.outdir,
            prefix=args.prefix, outgroups=args.outgroup)
    if args.command == "place":
        return cansnp_tree.place_main(
            scheme_dir=args.scheme_dir, organism=args.organism, tree_file=args.tree,
            alignment=args.alignment, reference=args.reference, genome=args.genome,
            members=args.members, min_match=args.min_match, output=args.output,
            report=args.report)
    if args.command == "trf":
        return cansnp_trf.run_trf(
            args.reference, args.output, args.match, args.mismatch, args.delta, args.pm,
            args.pi, args.min_score, args.max_period)
    parser.print_help()
    return 2


if __name__ == "__main__":
    sys.exit(main())
