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
from cansnpscorer.scoring_methods import ScoringManager, NeighborDistanceScoring, GCContentScoring, TandemRepeatProximityScoring, GeneProximityScoring, MutationTypeScoring

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


    # Initialize data
    snps = snp_file_parser(snp_file)
    reference = reference_parser(reference_file)
    genes = gff_parser(gff_file)
    tandem_repeats = trf_parser(trf_file)
    output = output_file

    # Initialize ScoringManager and Scoring methods
    scoring_methods = [NeighborDistanceScoring(), GCContentScoring(), TandemRepeatProximityScoring(), GeneProximityScoring(), MutationTypeScoring()]
    scoring_manager = ScoringManager(scoring_methods)

    # Calculate and sort scores
    scores = scoring_manager.calculate_scores(snps, reference, genes, tandem_repeats)
    sorted_scores = sorted(scores.items(), key=itemgetter(1), reverse=True)

    # Write output
    score_file_writer(sorted_scores, output)


if __name__ == "__main__":
    main()
