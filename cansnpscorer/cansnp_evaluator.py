import sys
import csv
import os
import argparse
from operator import itemgetter
from datetime import datetime
import logging
import shutil


from cansnpscorer.helpers.data_utils import filter_scores_by_cutoff, summarize_scores
from cansnpscorer.helpers.file_utils import reference_parser, gff_parser, trf_parser, snp_file_parser, score_file_writer, read_yaml_file, extract_yaml_data
from cansnpscorer.helpers.scoring_methods import NeighborDistanceScoring, GCContentScoring, TandemRepeatProximityScoring, GeneProximityScoring, MutationTypeScoring, ScoringManager





def main(snp_file=None, reference_file=None, gff_file=None, trf_file=None, hcfile=None, output_file=None, tandem_limit=None, gene_limit=None, gc_goal=None, cutoff_value=None, top_n_snp=None):
    if snp_file is None and reference_file is None:  # If no arguments passed in, use argparse
        parser = argparse.ArgumentParser(description="Evaluate CanSNP data with more contextual information.")
        parser.add_argument("--snp-file", type=str, required=True, help="Path to the input SNP file.")
        parser.add_argument("--reference-file", type=str, required=True, help="Path to the input reference genome file.")
        parser.add_argument("--gff-file", type=str, required=True, help="Path to the input GFF (General Feature Format) file.")
        parser.add_argument("--trf-file", type=str, required=True, help="Path to the input TRF (Tandem Repeats Finder) file.")
        parser.add_argument("--hcfile", type=str, required=True, help="Path to the input Hierarchical Clustering (HC) file.")
        parser.add_argument("--output-file", type=str, required=True, help="Path to the output file where results will be saved.")
        parser.add_argument("--tandem-limit", type=int, default=50, help="Tandem repeat proximity limit.")
        parser.add_argument("--gene-limit", type=int, default=0, help="Gene proximity limit.")
        parser.add_argument("--gc-goal", type=float, default=0.5, help="Goal for GC content scoring.")
        parser.add_argument("--cutoff-value", type=int, default=500, help="Cutoff value for SNP score filtering.")
        parser.add_argument("--top-n-snp", type=int, default=None, help="The number of top scoring SNPs to select for each node_id.")

        args = parser.parse_args()
        snps = snp_file_parser(args.snp_file)
        reference = reference_parser(args.reference_file)
        ref_gb = gff_parser(args.gff_file)
        trf = trf_parser(args.trf_file)
        output = args.output_file
        tandem_limit = args.tandem_limit
        gene_limit = args.gene_limit
        gc_goal = args.gc_goal
        cutoff_value = args.cutoff_value
        top_n_snp = args.top_n_snp

    # Initialize logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    snps = snp_file_parser(snp_file)
    reference = reference_parser(reference_file)
    ref_gb = gff_parser(gff_file)
    trf = trf_parser(trf_file)
    output = output_file

    output = output_file
    if os.path.isfile(output):
        logging.warning(f"Output file {output} already exists. Creating a backup.")
        timestamp = datetime.now().strftime('%Y%m%d%H%M%S')
        backup_name = f"{output}_{timestamp}.backup"
        shutil.copy(output, backup_name)
        logging.info(f"Backup created: {backup_name}")


    # Get the values from argparse


    correction = 0
    chromosome = 1
    chromosome_length = reference[1]

    scores = []
    for i in range(0, len(snps)):
        if snps[i][1] == "None":
            continue

        if snps[i][0] - correction > chromosome_length:
            correction = reference[chromosome]
            chromosome += 1

        final_score = ScoringManager(i, snps, reference, ref_gb, trf, chromosome_length, correction, chromosome,
                                     tandem_limit, gene_limit, gc_goal, ['gene_score'])

        scores.append(final_score)


    tree_data = read_yaml_file(hcfile)
    nodes_data = extract_yaml_data(tree_data)

    # Sort and write scores to file
    scores = sorted(scores, key=lambda x: x[0], reverse=False)
    filtered_scores = filter_scores_by_cutoff(scores, cutoff_value, top_n_snp)

    summarize_scores(nodes_data, filtered_scores)


    score_file_writer(filtered_scores, output)


if __name__ == "__main__":
    main()
