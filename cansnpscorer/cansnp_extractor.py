import argparse
import logging
from datetime import datetime
import os


from cansnpscorer.helpers.file_utils import fasta_parser, HCParser, cansnp_file_writer, read_yaml_file, recognize_hc_file_format, convert_data, save_to_file, recognize_group_file_format, read_hc_file
from cansnpscorer.helpers.data_utils import snp_extractor, canSNP_extractor


def main(groupfile=None, hcfile=None, fastafile=None, output=None):

    # Initialize logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    logging.info(f"Started analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info(f"Parameters: groupfile={groupfile}, hcfile={hcfile}, fastafile={fastafile}, output={output}")

    # Check if the input files exist
    if not (os.path.isfile(groupfile) and os.path.isfile(hcfile) and os.path.isfile(fastafile)):
        logging.error(f"Input files do not exist. Exiting...")
        sys.exit(1)

    # Check the format of the hc file
    hc_format = recognize_hc_file_format(hcfile)

    tree_data = read_hc_file(hcfile, hc_format)

    # Parse hierarchical data from file
    group_format = recognize_group_file_format(groupfile)

    hc_parser_obj = HCParser(groupfile, data_format=group_format)
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
