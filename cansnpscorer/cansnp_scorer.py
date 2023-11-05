import argparse
import logging
from datetime import datetime

from cansnpscorer import cansnp_evaluator, cansnp_extractor, cansnp_trf, cansnp_converter

# Import the version number
from cansnpscorer import __version__

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
    parser_extract.add_argument("--fastafile", required=True, help="Path to the input FASTA file containing sequence alignments.")
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
    parser_evaluate.add_argument("--cutoff-value", type=int, default=500, help="Cutoff value for SNP score filtering.")
    parser_evaluate.add_argument("--top-n-snp", type=int, default=None, help="The number of top scoring SNPs to select for each node_id.")
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

    # Create the parser for the "convert" command
    parser_convert = subparsers.add_parser('convert', help='Convert hierarchical data formats')
    parser_convert.add_argument("-i", "--infile", type=str, help="The input file name")
    parser_convert.add_argument("-p", "--prefix", type=str, help="The output file name")
    parser_convert.add_argument("--format_in", type=str, default=None, help="The format of the input file ('yaml', 'json', 'newick', 'child_parent')")
    parser_convert.add_argument("--format_out", type=str, required=True, help="The format of the output file ('yaml', 'json', 'newick', 'child_parent')")
    parser_convert.add_argument('--version', action='version', version=f'%(prog)s {__version__}', help="Show the current version of the script.")


    # Parse the arguments
    args = parser.parse_args()

    # Initialize logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    if args.command == 'extract':
        logging.info(f"Starting canSNP extraction...")
        cansnp_extractor.main(args.groupfile, args.hcfile, args.fastafile, args.output)
    elif args.command == 'evaluate':
        logging.info(f"Starting canSNP evaluation...")
        cansnp_evaluator.main(args.snp_file, args.reference_file, args.gff_file, args.trf_file, args.hcfile, args.output_file, args.tandem_limit, args.gene_limit, args.gc_goal, args.cutoff_value, args.top_n_snp)
    elif args.command == 'trf':
        logging.info(f"Starting canSNP TRF...")
        cansnp_trf.run_trf(
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
    elif args.command == 'convert':
        logging.info(f"Starting hierarchical data conversion...")
        cansnp_converter.main(args.infile, args.prefix, args.format_in, args.format_out)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
