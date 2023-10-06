import sys
import csv
import os
import argparse
from operator import itemgetter
from datetime import datetime
import logging
import shutil


from cansnpscorer.helpers.data_utils import neighbour_distance, gc_content, search, gene_checker
from cansnpscorer.helpers.file_utils import reference_parser, gff_parser, trf_parser, snp_file_parser, score_file_writer

def main(snp_file=None, reference_file=None, gff_file=None, trf_file=None, output_file=None):
    if snp_file is None and reference_file is None:  # If no arguments passed in, use argparse
        parser = argparse.ArgumentParser(description="Evaluate CanSNP data with more contextual information.")
        parser.add_argument("--snp-file", type=str, required=True, help="Path to the input SNP file.")
        parser.add_argument("--reference-file", type=str, required=True, help="Path to the input reference genome file.")
        parser.add_argument("--gff-file", type=str, required=True, help="Path to the input GFF (General Feature Format) file.")
        parser.add_argument("--trf-file", type=str, required=True, help="Path to the input TRF (Tandem Repeats Finder) file.")
        parser.add_argument("--output-file", type=str, required=True,
                            help="Path to the output file where results will be saved.")
        args = parser.parse_args()

        snp_file = args.snp_file
        reference_file = args.reference_file
        gff_file = args.gff_file
        trf_file = args.trf_file
        output_file = args.output_file


    # Initialize logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # snp_file = args.snp_file
    # reference_file = args.reference_file
    # gff_file = args.gff_file
    # trf_file = args.trf_file
    # output_file = args.output_file


    snps = snp_file_parser(snp_file)
    reference = reference_parser(reference_file)
    ref_gb = gff_parser(gff_file)
    trf = trf_parser(trf_file)
    output = output_file

    if os.path.isfile(output):
        logging.warning(f"Output file {output} already exists. Creating a backup.")
        timestamp = datetime.now().strftime('%Y%m%d%H%M%S')
        backup_name = f"{output}_{timestamp}.backup"
        shutil.copy(output, backup_name)
        logging.info(f"Backup created: {backup_name}")


    modifier = 100

    tandem_limit = 50
    gene_limit = 0
    gc_goal = 0.5

    correction = 0
    chromosome = 1

    chromosome_length = reference[1]
    scores = []
    for i in range(0,len(snps)):
        if snps[i][1] == "None":
            continue

        if snps[i][0] - correction > chromosome_length:
            correction = reference[chromosome]
            chromosome += 1
        if i > 0 and i < len(snps) - 1:
            neighbour = neighbour_distance(snps[i][0], snps[i-1][0], snps[i+1][0])

        elif i == 0:
            neighbour = neighbour_distance(snps[i][0], 0, snps[i+1])

        elif i == len(snps):
            neighbour = neighbour_distance(snps[i][0], snps[i-1], len(reference[0]))
        gc = gc_content(snps[i][0], reference[0], gc_goal)
        snp_pos = snps[i][0] - correction
        tandem = search(snp_pos, trf[chromosome], tandem_limit)
        gene = search(snp_pos, ref_gb[chromosome], gene_limit)
        mutation = 0

        if gene[0] == 0:
            mutation = gene_checker(snps[i],reference, gene[1], correction)
        score = [neighbour*1, gc*1 ,tandem[0]*1, (1 - gene[0])*2, mutation*1]
        score = [int(x*modifier) for x in score]
        total =  sum(score)

        scores.append([snps[i][0],snps[i][1], total] + score)

    scores = sorted(scores, key=itemgetter(2), reverse = True)
    score_file_writer(scores,output)

if __name__ == "__main__":
    main()
