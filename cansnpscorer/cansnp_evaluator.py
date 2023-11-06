import sys
import csv
import os
import argparse

from datetime import datetime
import logging
import shutil


from cansnpscorer.helpers.data_utils import filter_scores_by_cutoff, summarize_scores
from cansnpscorer.helpers.file_utils import reference_parser, gff_parser, trf_parser, snp_file_parser, score_file_writer, read_yaml_file, extract_yaml_data
from cansnpscorer.helpers.scoring_methods import NeighborDistanceScoring, GCContentScoring, TandemRepeatProximityScoring, GeneProximityScoring, MutationTypeScoring, ScoringManager





def main(snp_file=None, reference_file=None, gff_file=None, trf_file=None, hcfile=None, output_file=None, tandem_limit=None, gene_limit=None, gc_goal=None, cutoff_value=None, top_n_snp=None):

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

