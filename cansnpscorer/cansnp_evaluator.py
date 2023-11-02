import sys
import csv
import os
import argparse
from operator import itemgetter
from datetime import datetime
import logging
import shutil



from cansnpscorer.helpers.file_utils import reference_parser, gff_parser, trf_parser, snp_file_parser, score_file_writer
from cansnpscorer.helpers.scoring_methods import NeighborDistanceScoring, GCContentScoring, TandemRepeatProximityScoring, GeneProximityScoring, MutationTypeScoring, ScoringManager


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


    tandem_limit = 50
    gene_limit = 0
    gc_goal = 0.5

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

    # Sort and write scores to file
    scores = sorted(scores, key=lambda x: x[0], reverse=False)
    score_file_writer(scores, output)


if __name__ == "__main__":
    main()
