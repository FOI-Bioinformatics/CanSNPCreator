import logging
import os
from datetime import datetime

from cansnpscorer.helpers.data_utils import canSNP_extractor, snp_extractor
from cansnpscorer.helpers.file_utils import (
    HCParser,
    cansnp_file_writer,
    read_alignment,
    read_hc_file,
    recognize_group_file_format,
    recognize_hc_file_format,
)


def main(groupfile=None, hcfile=None, fastafile=None, output=None) -> int:
    """Extract canSNPs from an alignment. Returns a process exit status."""
    logging.info(f"Started analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info(f"Parameters: groupfile={groupfile}, hcfile={hcfile}, fastafile={fastafile}, output={output}")

    missing = [path for path in (groupfile, hcfile, fastafile) if not os.path.isfile(path)]
    if missing:
        logging.error(f"These input files do not exist: {', '.join(missing)}")
        return 1

    # Check the format of the hc file
    hc_format = recognize_hc_file_format(hcfile)

    tree_data = read_hc_file(hcfile, hc_format)

    # Parse hierarchical data from file
    group_format = recognize_group_file_format(groupfile)

    hc_parser_obj = HCParser(groupfile, data_format=group_format)
    hc_data = hc_parser_obj.reader(tree_data)


    # Parse the alignment, which must hold the reference as its first record
    try:
        alignment = read_alignment(fastafile)
    except (OSError, ValueError) as error:
        logging.error(error)
        return 1

    # Extract SNPs
    snps = snp_extractor(alignment)

    # ExtractCanSNPs
    cansnps = canSNP_extractor(alignment, hc_data, snps)

    # Write canSNPs to the output file
    cansnp_file_writer(cansnps, output)

    logging.info(f"Finished analysis at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    return 0
