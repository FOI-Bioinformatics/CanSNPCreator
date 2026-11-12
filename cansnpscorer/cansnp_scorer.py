import argparse
import logging
import sys

from cansnpscorer import (
    __version__,
    cansnp_audit,
    cansnp_evaluator,
    cansnp_export,
    cansnp_extractor,
    cansnp_tree,
    cansnp_trf,
)
from cansnpscorer.helpers.data_utils import MINIMUM_SNPS_PER_NODE
from cansnpscorer.helpers.scoring_methods import DEFAULT_WEIGHTS, MAXIMUM_SCORE, TERMS


def score_weights(settings):
    """Turn --score-weight TERM=FACTOR arguments into a dictionary of weights."""
    if not settings:
        return dict(DEFAULT_WEIGHTS)

    weights = {}
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


def main():
    # Create the top-level parser
    parser = argparse.ArgumentParser(
        prog="CanSNPCreator",
        description="Extract, evaluate and export canonical SNPs from a whole genome alignment.")
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}',
                                help="Show the current version of the script.")
    parser.add_argument('--log-level', default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        help="How much to report while running. Default: INFO.")
    subparsers = parser.add_subparsers(dest='command', help='Sub-command help')

    # Create the parser for the "extract" command
    parser_extract = subparsers.add_parser('extract', help='Extract canSNPs')
    parser_extract.add_argument("--groupfile", required=True, help="Path to the input group file containing genome data.")
    parser_extract.add_argument("--hcfile", required=True, help="Path to the input Hierarchical Clustering (HC) file.")
    parser_extract.add_argument("--fastafile", required=True, help="Path to the input FASTA file containing sequence alignments. Its first record must be the reference genome.")
    parser_extract.add_argument("--reference-file", required=True,
                                help="Path to the reference genome in FASTA format, one record "
                                     "per contig, in the order the alignment concatenates them. "
                                     "It gives the contig names and boundaries that positions "
                                     "are reported against.")
    parser_extract.add_argument("-o", "--output", required=True, help="Path to the output file where extracted SNPs and canSNPs will be saved.")
    parser_extract.add_argument('--version', action='version', version=f'%(prog)s {__version__}', help="Show the current version of the script.")

    # Create the parser for the "evaluate" command
    parser_evaluate = subparsers.add_parser('evaluate', help="Evaluate CanSNP data with more contextual information.")
    parser_evaluate.add_argument("--snp-file", type=str, required=True, help="Path to the input SNP file.")
    parser_evaluate.add_argument("--reference-file", type=str, required=True, help="Path to the input reference genome file.")
    parser_evaluate.add_argument("--gff-file", type=str, required=True, help="Path to the input GFF (General Feature Format) file.")
    parser_evaluate.add_argument("--trf-file", type=str, required=True, help="Path to the input TRF (Tandem Repeats Finder) file.")
    parser_evaluate.add_argument("--hcfile", type=str, required=True, help="Path to the input Hierarchical Clustering (HC) file.")
    parser_evaluate.add_argument("--output-file", type=str, required=True, help="Path to the output file where results will be saved.")
    parser_evaluate.add_argument("--tandem-limit", type=int, default=50, help="Tandem repeat proximity limit.")
    parser_evaluate.add_argument("--gene-limit", type=int, default=0, help="Gene proximity limit.")
    parser_evaluate.add_argument("--gc-goal", type=float, default=0.5, help="Goal for GC content scoring.")
    parser_evaluate.add_argument("--cutoff-value", type=int, default=500,
                                 help=f"Only report SNPs scoring above this value. The highest "
                                      f"reachable total is {MAXIMUM_SCORE} with the default "
                                      f"weights, and the default cutoff of 500 in practice "
                                      f"requires the SNP to fall inside a gene.")
    parser_evaluate.add_argument("--top-n-snp", type=int, default=None, help="The number of top scoring SNPs to select for each node_id.")
    parser_evaluate.add_argument("--score-weight", action="append", metavar="TERM=FACTOR",
                                 help=f"Weight one score term before the terms are summed, for "
                                      f"example gene_score=2. Terms: {', '.join(TERMS)}. "
                                      f"Repeatable. Default: "
                                      f"{', '.join(f'{k}={v:g}' for k, v in DEFAULT_WEIGHTS.items())}.")
    parser_evaluate.add_argument('--version', action='version', version=f'%(prog)s {__version__}', help="Show the current version of the script.")

    # Create the parser for the "export" command
    parser_export = subparsers.add_parser(
        'export', help="Write a scored scheme as CanSNPer2-data source files.")
    parser_export.add_argument("--score-file", required=True,
                               help="The CSV file the evaluate subcommand wrote.")
    parser_export.add_argument("--hcfile", required=True,
                               help="Path to the input Hierarchical Clustering (HC) file.")
    parser_export.add_argument("--reference-file", required=True,
                               help="The reference genome, whose contig lengths turn positions "
                                    "into the concatenated coordinates the SNP file uses.")
    parser_export.add_argument("--organism", required=True,
                               help="Names the three files, as in francisella_tularensis.")
    parser_export.add_argument("--genome", required=True,
                               help="The short name of the reference genome, as the SNP file "
                                    "refers to it. It may not hold an underscore.")
    parser_export.add_argument("--strain", required=True, help="The strain column of both files.")
    parser_export.add_argument("--reference", dest="reference_label", required=True,
                               help="What the SNPs are cited as, as in Svensson2009.")
    parser_export.add_argument("--genbank-id", required=True, help="As in GCA_000007125.1.")
    parser_export.add_argument("--refseq-id", required=True, help="As in GCF_000007125.1.")
    parser_export.add_argument("--assembly-name", required=True, help="As in ASM712v1.")
    parser_export.add_argument("--outdir", default=".", help="Where to write the files.")
    parser_export.add_argument("--prune", action="store_true",
                               help="Leave out nodes that have no canSNP, instead of refusing to "
                                    "write a scheme. Their children then hang from the nearest "
                                    "node that has one.")
    parser_export.add_argument('--version', action='version', version=f'%(prog)s {__version__}')

    # Create the parser for the "audit" command
    parser_audit = subparsers.add_parser(
        'audit', help="Examine an existing CanSNPer2-data scheme and report whether it should "
                      "be extended or modified.")
    parser_audit.add_argument("--scheme-dir", required=True,
                              help="The directory holding a scheme's three source files, as in "
                                   "the source directory of CanSNPer2-data.")
    parser_audit.add_argument("--organism", default=None,
                              help="Which scheme to read, as in francisella_tularensis. It may "
                                   "be left out when the directory holds one scheme.")
    parser_audit.add_argument("--reference-file", default=None,
                              help="The reference genome the scheme is annotated against, in "
                                   "FASTA format. With it, the base at every position is checked "
                                   "against the two alleles, which is the direct test for "
                                   "positions taken against another version of the assembly.")
    parser_audit.add_argument("--score-file", default=None,
                              help="A score file from the evaluate subcommand, to report where a "
                                   "better scoring canSNP is available. Needs --reference-file.")
    parser_audit.add_argument("--genome", default=None,
                              help="Which of the scheme's reference genomes the reference file "
                                   "holds. Needed only when the scheme names more than one.")
    parser_audit.add_argument("--min-required-hits", type=int, default=MINIMUM_SNPS_PER_NODE,
                              help=f"The number of SNPs CanSNPer2 requires on the path to a node "
                                   f"before it calls it. Default: {MINIMUM_SNPS_PER_NODE}, as "
                                   f"CanSNPer2 itself uses.")
    parser_audit.add_argument("--proximity", type=int, default=50,
                              help="Report two canSNPs closer than this many bases, since a "
                                   "neighbouring SNP falls inside an assay. Default: 50.")
    parser_audit.add_argument("--fastafile", default=None,
                              help="A whole genome alignment with the reference genome as its "
                                   "first record, either gapped or projected onto reference "
                                   "coordinates as RepGenR2 and snippy write it. With "
                                   "--groupfile, each SNP is checked in every genome, and the "
                                   "genomes are searched for clades the scheme has no node for. "
                                   "Needs --reference-file.")
    parser_audit.add_argument("--groupfile", default=None,
                              help="Where the alignment's genomes sit in the scheme: a Strain "
                                   "and Group CSV, a repgenr group file, or two tab separated "
                                   "columns of genome and node, as RepGenR2's genomes_map.tsv. "
                                   "A genome may be placed at an internal node.")
    parser_audit.add_argument("--tree", default=None,
                              help="Instead of --groupfile: a rooted genome tree whose leaves are "
                                   "the alignment's genomes. The genomes are placed on the "
                                   "scheme by the tree, as the place subcommand does, rather "
                                   "than by typing them with the SNPs under test.")
    parser_audit.add_argument("--genome-map", default=None,
                              help="With --tree: RepGenR2's genomes_map.tsv.")
    parser_audit.add_argument("--min-clade-size", type=int, default=2,
                              help="Report a clade the scheme does not resolve only when it "
                                   "holds at least this many genomes. Default: 2.")
    parser_audit.add_argument("--min-clade-snps", type=int, default=1,
                              help="Report a clade the scheme does not resolve only when at "
                                   "least this many SNPs mark it. Default: 1.")
    parser_audit.add_argument("-o", "--output", default=None,
                              help="Write the findings to this tab separated file as well.")
    parser_audit.add_argument("--strict", action="store_true",
                              help="Exit 2 for a recommendation as well as for a problem, for "
                                   "use in continuous integration.")
    parser_audit.add_argument('--version', action='version', version=f'%(prog)s {__version__}')

    # Create the parser for the "collapse" command
    parser_collapse = subparsers.add_parser(
        'collapse', help="Reduce a genome tree to the nodes of a scheme, and place its genomes.")
    parser_collapse.add_argument("--tree", required=True,
                                 help="A rooted genome tree whose leaves are genomes: Newick, "
                                      "or a child and parent table such as RepGenR2's "
                                      "tree2tax.tsv.")
    parser_collapse.add_argument("--nodes-file", default=None,
                                 help="RepGenR2's tree2tax_nodes.tsv, for the support and "
                                      "branch length of each node. A Newick tree carries both.")
    parser_collapse.add_argument("--genome-map", default=None,
                                 help="RepGenR2's genomes_map.tsv, so that genomes folded into a "
                                      "representative by dereplication are placed with it.")
    parser_collapse.add_argument("--min-support", type=float, default=None,
                                 help="Remove internal nodes with a lower support, as a fraction; "
                                      "percentages in the tree are read as fractions.")
    parser_collapse.add_argument("--min-genomes", type=int, default=2,
                                 help="Remove internal nodes with fewer genomes below them. "
                                      "Default: 2.")
    parser_collapse.add_argument("--min-branch-length", type=float, default=None,
                                 help="Remove internal nodes whose branch is shorter, in the "
                                      "tree's units.")
    parser_collapse.add_argument("--outgroup", action="append", default=[], metavar="LEAF",
                                 help="A leaf of the tree that is an outgroup: it is placed at "
                                      "no node, and stays in the alignment as a genome outside "
                                      "every node, which gives the root a polarity. Repeatable.")
    parser_collapse.add_argument("--rename", default=None, metavar="PREFIX",
                                 help="Name the nodes PREFIX, PREFIX.1, PREFIX.1.2 and so on, "
                                      "largest first, instead of the tree's own names.")
    parser_collapse.add_argument("--outdir", default=".", help="Where to write the files.")
    parser_collapse.add_argument("--prefix", default="collapsed",
                                 help="Names the three files: PREFIX_hierarchy.tsv, "
                                      "PREFIX_genomes.tsv and PREFIX_nodes.tsv.")
    parser_collapse.add_argument('--version', action='version', version=f'%(prog)s {__version__}')

    # Create the parser for the "place" command
    parser_place = subparsers.add_parser(
        'place', help="Place the genomes of a tree on the nodes of an existing scheme.")
    parser_place.add_argument("--scheme-dir", required=True,
                              help="The directory holding the scheme's three source files.")
    parser_place.add_argument("--organism", default=None, help="Which scheme to read.")
    parser_place.add_argument("--tree", required=True,
                              help="A rooted genome tree whose leaves are the alignment's "
                                   "genomes, preferably rooted on an outgroup.")
    parser_place.add_argument("--fastafile", required=True,
                              help="The whole genome alignment, reference first.")
    parser_place.add_argument("--reference-file", required=True,
                              help="The reference genome the scheme is annotated against.")
    parser_place.add_argument("--genome", default=None,
                              help="Which of the scheme's reference genomes that is, when it "
                                   "names more than one.")
    parser_place.add_argument("--genome-map", default=None,
                              help="RepGenR2's genomes_map.tsv, so that genomes folded into a "
                                   "leaf are placed with it.")
    parser_place.add_argument("--min-match", type=float, default=0.5,
                              help="The lowest Jaccard index between a node's derived genomes "
                                   "and a clade that counts as a match. Default: 0.5.")
    parser_place.add_argument("-o", "--output", required=True,
                              help="The genome map to write: genome and node, tab separated.")
    parser_place.add_argument("--report", default=None,
                              help="Also write how each node was matched to the tree.")
    parser_place.add_argument('--version', action='version', version=f'%(prog)s {__version__}')

    # Create the parser for the "trf" command
    parser_trf = subparsers.add_parser('trf', help='Perform TRF on canSNPs')
    parser_trf.add_argument('--input', required=True, help='Input file for TRF')
    parser_trf.add_argument('--output', required=True, help='Output file for TRF results')
    parser_trf.add_argument('--match', type=int, default=2, help='Matching weight. Default is 2.')
    parser_trf.add_argument('--mismatch', type=int, default=7, help='Mismatching penalty. Default is 7.')
    parser_trf.add_argument('--delta', type=int, default=7, help='Indel penalty. Default is 7.')
    parser_trf.add_argument('--pm', type=int, default=80, help='Match probability (whole number). Default is 80.')
    parser_trf.add_argument('--pi', type=int, default=10, help='Indel probability (whole number). Default is 10.')
    parser_trf.add_argument('--minscore', type=int, default=10, help='Minimum alignment score to report. Default is 10.')
    parser_trf.add_argument('--max_period', type=int, default=2000, help='Maximum period size to report. Default is 2000.')
    parser_trf.add_argument('--version', action='version', version=f'%(prog)s {__version__}', help="Show the current version of the script.")

    # Parse the arguments
    args = parser.parse_args()

    # Initialize logging
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format='%(asctime)s - %(levelname)s - %(message)s')

    if args.command == 'extract':
        logging.info("Starting canSNP extraction")
        return cansnp_extractor.main(args.groupfile, args.hcfile, args.fastafile,
                                     args.reference_file, args.output)
    elif args.command == 'evaluate':
        logging.info("Starting canSNP evaluation")
        try:
            weights = score_weights(args.score_weight)
        except argparse.ArgumentTypeError as error:
            parser_evaluate.error(str(error))
        return cansnp_evaluator.main(
            args.snp_file, args.reference_file, args.gff_file, args.trf_file, args.hcfile,
            args.output_file, args.tandem_limit, args.gene_limit, args.gc_goal,
            args.cutoff_value, args.top_n_snp, weights)
    elif args.command == 'export':
        logging.info("Writing the scheme")
        return cansnp_export.main(
            args.score_file, args.hcfile, args.reference_file, args.outdir, args.organism,
            args.genome, args.strain, args.reference_label, args.genbank_id, args.refseq_id,
            args.assembly_name, args.prune)
    elif args.command == 'audit':
        logging.info("Auditing the scheme")
        return cansnp_audit.main(
            args.scheme_dir, args.organism, args.reference_file, args.score_file, args.genome,
            args.min_required_hits, args.proximity, args.output, args.strict,
            args.fastafile, args.groupfile, args.min_clade_size, args.min_clade_snps,
            args.tree, args.genome_map)
    elif args.command == 'collapse':
        logging.info("Collapsing the genome tree")
        return cansnp_tree.collapse_main(
            args.tree, args.nodes_file, args.genome_map, args.min_support, args.min_genomes,
            args.min_branch_length, args.rename, args.outdir, args.prefix, args.outgroup)
    elif args.command == 'place':
        logging.info("Placing the genomes on the scheme")
        return cansnp_tree.place_main(
            args.scheme_dir, args.organism, args.tree, args.fastafile, args.reference_file,
            args.genome, args.genome_map, args.min_match, args.output, args.report)
    elif args.command == 'trf':
        logging.info("Starting canSNP TRF")
        return cansnp_trf.run_trf(
            args.input,
            args.output,
            args.match,
            args.mismatch,
            args.delta,
            args.pm,
            args.pi,
            args.minscore,
            args.max_period
        )
    else:
        parser.print_help()
        return 2


if __name__ == "__main__":
    sys.exit(main())
