import argparse
import logging
from datetime import datetime
import os


from cansnpscorer.helpers.file_utils import fasta_parser, HCParser, cansnp_file_writer, read_yaml_file, recognize_hc_file_format, convert_data, save_to_file, recognize_group_file_format, read_hc_file
from cansnpscorer.helpers.data_utils import snp_extractor, canSNP_extractor


def main(csvfile=None, hcfile=None, fastafile=None, output=None):
    if csvfile is None and hcfile is None:  # If no arguments passed in, use argparse
        parser = argparse.ArgumentParser(
            description="Extract Single Nucleotide Polymorphisms (SNPs) and Canonical SNPs (canSNPs) from genomic data.")
        parser.add_argument("--csvfile", required=True, help="Path to the input CSV file containing genome data.")
        parser.add_argument("--hcfile", required=True, help="Path to the input Hierarchical Clustering (HC) file.")
        parser.add_argument("--fastafile", required=True,
                            help="Path to the input FASTA file containing sequence alignments.")
        parser.add_argument("-o", "--output", required=True,
                            help="Path to the output file where extracted SNPs and canSNPs will be saved.")
        args = parser.parse_args()

        csvfile = args.csvfile
        hcfile = args.hcfile
        fastafile = args.fastafile
        output = args.output

    # Initialize logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    logging.info(f"Started analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info(f"Parameters: csvfile={csvfile}, hcfile={hcfile}, fastafile={fastafile}, output={output}")

    # Check if the input files exist
    if not (os.path.isfile(csvfile) and os.path.isfile(hcfile) and os.path.isfile(fastafile)):
        logging.error(f"Input files do not exist. Exiting...")
        sys.exit(1)

    # Check the format of the hc file
    hc_format = recognize_hc_file_format(hcfile)

    tree_data = read_hc_file(hcfile, hc_format)

    # Parse hierarchical data from file
    group_format = recognize_group_file_format(csvfile)

    hc_parser_obj = HCParser(csvfile, data_format=group_format)
    hc_data = hc_parser_obj.reader(tree_data)

    # Parse FASTA data
    fasta_data = fasta_parser(fastafile, 'dict')

    # Extract SNPs
    snps = snp_extractor(fasta_data)

    # ExtractCanSNPs
    cansnps = canSNP_extractor(fasta_data, hc_data, snps)

    # Write canSNPs to the output file
    cansnp_file_writer(cansnps, output)

    logging.info(f"Finished analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()
