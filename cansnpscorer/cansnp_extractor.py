import argparse
import yaml

from cansnpscorer.helpers.file_utils import fasta_parser, hc_parser, cansnp_file_writer
from cansnpscorer.helpers.data_utils import snp_extractor, canSNP_extractor

def main():
    parser = argparse.ArgumentParser(
        description="Extract Single Nucleotide Polymorphisms (SNPs) and Canonical SNPs (canSNPs) from genomic data.")

    parser.add_argument("--csvfile", required=True, help="Path to the input CSV file containing genome data.")
    parser.add_argument("--hcfile", required=True, help="Path to the input Hierarchical Clustering (HC) file.")
    parser.add_argument("--fastafile", required=True,
                        help="Path to the input FASTA file containing sequence alignments.")
    parser.add_argument("-o", "--output", required=True,
                        help="Path to the output file where extracted SNPs and canSNPs will be saved.")

    args = parser.parse_args()

    # Parse HC data from YAML file
    hc_parser_obj = hc_parser(args.csvfile)

    with open(args.hcfile, 'r') as stream:
        yaml_data = next(yaml.load_all(stream, Loader=yaml.SafeLoader))

    hc_data = hc_parser_obj.reader(yaml_data)

    # Parse FASTA data
    fasta_data = fasta_parser(args.fastafile, True, True)

    # Extract SNPs and canSNPs
    snps = snp_extractor(fasta_data)
    cansnps = canSNP_extractor(fasta_data, hc_data, snps)

    # Write canSNPs to the output file
    cansnp_file_writer(cansnps, args.output)

if __name__ == "__main__":
    main()
