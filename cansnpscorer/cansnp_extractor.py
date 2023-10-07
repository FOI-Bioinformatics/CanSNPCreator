import argparse
import logging
from datetime import datetime


from cansnpscorer.helpers.file_utils import fasta_parser, HCParser, cansnp_file_writer, read_yaml_file
from cansnpscorer.helpers.data_utils import snp_extractor, canSNP_extractor


def main(csvfile=None, hcfile=None, fastafile=None, output=None, num_threads=4):
    if csvfile is None and hcfile is None:  # If no arguments passed in, use argparse
        parser = argparse.ArgumentParser(
            description="Extract Single Nucleotide Polymorphisms (SNPs) and Canonical SNPs (canSNPs) from genomic data.")
        parser.add_argument("--csvfile", required=True, help="Path to the input CSV file containing genome data.")
        parser.add_argument("--hcfile", required=True, help="Path to the input Hierarchical Clustering (HC) file.")
        parser.add_argument("--fastafile", required=True,
                            help="Path to the input FASTA file containing sequence alignments.")
        parser.add_argument("-o", "--output", required=True,
                            help="Path to the output file where extracted SNPs and canSNPs will be saved.")
        parser.add_argument("-t", "--threads", type=int, default=4,
                            help="Number of threads to use for parallel processing. Default is 4.")
        args = parser.parse_args()

        csvfile = args.csvfile
        hcfile = args.hcfile
        fastafile = args.fastafile
        output = args.output
        num_threads = args.threads

    # Initialize logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    logging.info(f"Started analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info(f"Parameters: csvfile={csvfile}, hcfile={hcfile}, fastafile={fastafile}, output={output}, threads={num_threads}")


    # Parse HC data from YAML file
    hc_parser_obj = HCParser(csvfile)

    yaml_data = read_yaml_file(hcfile)

    hc_data = hc_parser_obj.reader(yaml_data)

    # Parse FASTA data
    fasta_data = fasta_parser(fastafile, 'dict')

    # Extract SNPs
    snps = snp_extractor(fasta_data, num_threads)

    # ExtractCanSNPs
    cansnps = canSNP_extractor(fasta_data, hc_data, snps)

    # Write canSNPs to the output file
    cansnp_file_writer(cansnps, outputfile_name)

    logging.info(f"Finished analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()
