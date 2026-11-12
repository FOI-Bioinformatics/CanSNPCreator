import argparse
import logging
import sys

from cansnpscorer import __version__, cansnp_evaluator, cansnp_extractor, cansnp_trf
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
    parser = argparse.ArgumentParser(description='CanSNP Scorer')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}',
                                help="Show the current version of the script.")
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
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

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
