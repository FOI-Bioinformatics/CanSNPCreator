
import re
import os.path
import sys
import csv
import logging
from ruamel.yaml import YAML
from typing import List, Dict, Union, NoReturn, Tuple
from operator import itemgetter

def read_yaml_file(file: str) -> Union[Dict, None]:
    """
    Read and load the first YAML document from a file.

    Parameters:
        file (str): Path to the YAML file.

    Returns:
        Union[Dict, None]: Dictionary containing the loaded data, or None if an error occurs.
    """

    logging.info(f"Attempting to read YAML data from {file}")

    try:
        yaml = YAML(typ='safe')
        with open(file, 'r') as yam:
            yaml_data = next(yaml.load_all(yam))
        logging.info(f"Successfully read YAML data from {file}")
        return yaml_data
    except FileNotFoundError:
        logging.error(f"File {file} not found")
    except PermissionError:
        logging.error(f"Permission denied: Cannot read {file}")
    except Exception as e:
        logging.error(f"Failed to read YAML data from {file}. Exception: {e}")
    return None

def clean_name(name: str) -> str:
    """
    Cleans a given string by removing all characters that are not alphanumeric.

    Args:
        name (str): The string to clean.

    Returns:
        str: A cleaned string containing only alphanumeric characters.

    Example:
        >>> clean_name("Example String!")
        'ExampleString'
    """

    logging.debug(f"Cleaning the name: {name}")
    cleaned_name = re.sub('[^A-Za-z0-9]+', '', name)
    logging.debug(f"Cleaned name: {cleaned_name}")

    return cleaned_name



def fasta_parser(fasta_file: str, output_format: str = 'list') -> Union[List[str], List[Tuple[str, str]], Dict[str, str]]:
    """
    Parse a FASTA file.

    Args:
    fasta_file: Path to the FASTA file.
    output_format: Output format ('list', 'dict', 'named_list').

    Returns:
    Parsed genomes in the specified output format.
    """

    logging.info(f"Starting to parse the FASTA file {fasta_file} with output format {output_format}.")

    genomes = [] if output_format != 'dict' else {}
    processed_entries = 0

    try:
        with open(fasta_file, 'r') as f:
            name, genome = "", ""
            for line in f:
                if line.startswith(">"):
                    if name:
                        genome = genome.replace("\n", "")
                        if output_format == 'dict':
                            genomes[name] = genome
                        elif output_format == 'named_list':
                            genomes.append((name, genome))
                        else:
                            genomes.append(genome)
                        processed_entries += 1

                    name = clean_name(line[1:].strip())
                    genome = ""
                else:
                    genome += line

            if name:
                genome = genome.replace("\n", "")
                if output_format == 'dict':
                    genomes[name] = genome
                elif output_format == 'named_list':
                    genomes.append((name, genome))
                else:
                    genomes.append(genome)

            logging.info(f"Successfully parsed {processed_entries} entries from the FASTA file {fasta_file}.")
            return genomes

    except FileNotFoundError:
        logging.error(f"File {fasta_file} not found.")
        return None
    except PermissionError:
        logging.error(f"Permission denied. Cannot read {fasta_file}")
        return None
    except Exception as e:
        logging.error(f"An error occurred while parsing the FASTA file {fasta_file}. Exception: {e}")
        return None


def parse(fasta, names = False, dictonary = False):
	if not os.path.isfile(fasta):
		print("File input does not exist. Enter a correct filepath")
		return -1
	if isinstance(names,bool) and isinstance(dictonary, bool):
		return fasta_parser(fasta, names, dictonary)
	else:
		print("Second and third arguments must be booleans")
		return -1




class HCParser:
    """
    A class to handle hierarchical clustering data from CSV and YAML files.
    """
    def __init__(self, csvfile: str) -> None:
        """
        Initialize the HCParser object.

        Args:
            csvfile (str): The path to the CSV file containing group and strain information.
        """
        logging.info(f"Initializing HCParser with CSV file: {csvfile}")
        self.dictionary: Dict[str, List[str]] = {}
        self.group_label = "Group"
        self.group_index = 0
        self.strain_label = "Strain"
        self.strain_index = 0
        self.isolates: Dict[str, List[str]] = {}
        self._group_handler(csvfile)

    def _group_handler(self, csv_file: str) -> None:
        """
        Process the CSV file and populate the 'isolates' dictionary.

        Args:
            csv_file (str): The path to the CSV file.
        """
        logging.info(f"Processing CSV file: {csv_file}")
        with open(csv_file, 'r') as csvfile:
            reader = csv.reader(csvfile)
            for index, row in enumerate(reader):
                if index == 0:
                    if self.group_label in row and self.strain_label in row:
                        self.group_index = row.index(self.group_label)
                        self.strain_index = row.index(self.strain_label)
                    else:
                        logging.error("Group or Strain column not found in csv file. Exiting...")
                        sys.exit(1)
                else:
                    group = row[self.group_index]
                    strain = clean_name(row[self.strain_index])
                    self.isolates.setdefault(group, []).append(strain)

    def _dictionary_reader(self, yaml_data: Dict[str, Union[str, Dict]]) -> List[str]:
        """
        Recursively read YAML data and populate the 'dictionary'.

        Args:
            yaml_data (Dict): Dictionary containing YAML data.

        Returns:
            List[str]: A list of isolates.
        """
        isolates: List[str] = []
        name = list(yaml_data.keys())[0]
        items = list(yaml_data.values())[0]
        for group in items:
            if isinstance(group, dict):
                isolates += self._dictionary_reader(group)
            else:
                group_isolate = self.isolates.get(group, [])
                isolates += group_isolate
                self.dictionary[group] = group_isolate

        self.dictionary[name] = isolates
        return isolates

    def reader(self, lista: List[Union[str, Dict]]) -> Dict[str, List[str]]:
        """
        Read a list containing either strings or dictionaries and populate the 'dictionary'.

        Args:
            lista (List): The list containing either strings or dictionaries.

        Returns:
            Dict: The populated dictionary.
        """
        logging.info("Reading list to populate dictionary.")
        for item in lista:
            if isinstance(item, dict):
                self._dictionary_reader(item)
            else:
                self.dictionary[item] = self.isolates.get(item, [])

        return self.dictionary


def snp_file_parser(snpfile: str) -> List[List[Union[int, str]]]:
    """
    Parse a SNP file.

    Args:
    snpfile (str): The path to the SNP file.

    Returns:
    List[List[Union[int, str]]]: A list of lists, where each inner list represents a line in the SNP file
                                 and contains an integer representing the SNP position and a string for the variant.

    Raises:
    FileNotFoundError: If the SNP file is not found.
    PermissionError: If the permission to read the SNP file is denied.
    """

    logging.info(f"Starting to parse SNP file: {snpfile}")
    snps: List[List[Union[int, str]]] = []  # Explicitly mention the type for better readability

    try:
        with open(snpfile, 'r') as snp_file:
            for line_number, line in enumerate(snp_file, 1):  # Added line_number for better logging
                line = line.strip()
                if not line:
                    logging.debug(f"Ignoring empty line at line {line_number}.")
                    continue
                if line.startswith("#"):
                    logging.debug(f"Ignoring comment line at line {line_number}.")
                    continue

                variants = line.split(" ")
                try:
                    variants[0] = int(variants[0])
                    snps.append(variants)
                except ValueError:
                    logging.error(f"ValueError: Could not convert SNP {line} to an integer at line {line_number}. Ignoring SNP.")

        logging.info(f"Successfully parsed {len(snps)} SNPs from {snpfile}.")
        return snps
    except FileNotFoundError:
        logging.error(f"File {snpfile} not found.")
        return []
    except PermissionError:
        logging.error(f"Permission denied: Cannot read {snpfile}")
        return []

def reference_parser(reference: str) -> List[Union[str, int]]:
    """
    Parse a reference genome file.

    Args:
    reference (str): The path to the reference genome file.

    Returns:
    List[Union[str, int]]: A list containing the reference genome as a string and the lengths of each chromosome as integers.

    Raises:
    FileNotFoundError: If the reference file is not found.
    PermissionError: If the permission to read the reference file is denied.
    """

    logging.info(f"Starting to parse reference file: {reference}")
    ref: str = ""
    chromosomes: int = 0
    chro_lengths: List[int] = []

    try:
        with open(reference, 'r') as fil:
            for line_number, line in enumerate(fil, 1):
                line = line.strip()
                if not line:
                    logging.debug(f"Ignoring empty line at line {line_number}.")
                    continue

                if ">" not in line:
                    ref += line
                else:
                    if chromosomes != 0:
                        chro_lengths.append(len(ref))
                    chromosomes += 1
                    logging.debug(f"Found chromosome header at line {line_number}, total chromosomes now {chromosomes}.")

        return_list: List[Union[str, int]] = [ref] + chro_lengths
        logging.info(f"Successfully parsed {chromosomes} chromosomes from {reference}.")
        return return_list
    except FileNotFoundError:
        logging.error(f"File {reference} not found.")
        return []
    except PermissionError:
        logging.error(f"Permission denied: Cannot read {reference}")
        return []

def gff_parser(gff_file: str) -> Dict[int, List[List[Union[int, str]]]]:
    """
    Parse a GFF (General Feature Format) file to extract gene information.

    Args:
        gff_file (str): The path to the GFF file.

    Returns:
        Dict[int, List[List[Union[int, str]]]]: A dictionary where keys are region numbers and
                                                 values are lists of lists containing gene information.
                                                 Each inner list contains the start position, end position,
                                                 and type of the gene.

    Raises:
        FileNotFoundError: If the GFF file is not found.
        PermissionError: If the permission to read the GFF file is denied.
    """

    logging.info(f"Starting to parse GFF file: {gff_file}")
    genes: Dict[int, List[List[Union[int, str]]]] = {}

    try:
        with open(gff_file, 'r') as gff:
            current_region: int = 0
            for line_number, line in enumerate(gff, 1):
                line = line.strip()
                if "#" in line:
                    logging.debug(f"Ignoring comment line at line {line_number}.")
                    continue

                gene = line.split("\t")
                if "region" in gene:
                    current_region += 1
                    genes[current_region] = []
                    logging.debug(f"Found region header at line {line_number}, total regions now {current_region}.")
                elif "gene" in gene:
                    index = gene.index("gene")
                    genes[current_region].append([int(gene[index + 1]) - 1, int(gene[index + 2]) - 1, gene[index + 4]])

        logging.info(f"Successfully parsed {len(genes)} regions from {gff_file}.")
        return genes

    except FileNotFoundError:
        logging.error(f"File {gff_file} not found.")
        return {}
    except PermissionError:
        logging.error(f"Permission denied: Cannot read {gff_file}")
        return {}

def trf_parser(trf_file_path: str) -> Dict[int, List[Tuple[int, int]]]:
    """
    Parses a TRF (Tandem Repeats Finder) file and returns intervals for each chromosome.

    Args:
    trf_file_path: The file path of the TRF output.

    Returns:
    A dictionary where keys are chromosome numbers and values are lists of intervals.
    Each interval is represented as a tuple (start, end).
    """
    intervals_by_chromosome = {}
    chromosome_number = 0

    with open(trf_file_path, 'r') as trf_file:
        for line in trf_file:
            if "@" in line:
                chromosome_number += 1
                intervals_by_chromosome[chromosome_number] = []
                continue

            interval_data = line.strip("\n").split(" ")

            try:
                start, end = map(int, interval_data[:2])
                intervals_by_chromosome[chromosome_number].append((start, end))
            except ValueError:
                continue

    for chromosome in intervals_by_chromosome:
        intervals_by_chromosome[chromosome] = sorted(intervals_by_chromosome[chromosome], key=itemgetter(0))

    return intervals_by_chromosome


def score_file_writer(snp_data: List[List[Union[str, int, float]]], output_file_path: str) -> None:
    """
    Writes SNP data to a CSV file.

    Args:
    snp_data: A list of lists, where each inner list represents a row of SNP data.
    output_file_path: The file path where the CSV file will be saved.

    Returns:
    None
    """
    header = ["SNP", "Group", "Score", "Neighbour-snps", "GC-content", "Tandem", "Gene", "Mutation"]

    with open(output_file_path, 'w', newline='') as out_file:
        csv_writer = csv.writer(out_file, delimiter=',')
        csv_writer.writerow(header)
        csv_writer.writerows(snp_data)


def check_files_exist(files: List[str]) -> bool:
    """
    Check if all the specified files exist.

    Args:
        files (List[str]): A list of file paths to check for existence.

    Returns:
        bool: True if all files exist, otherwise False.

    Example:
        >>> check_files_exist(["file1.txt", "file2.txt"])
        True
    """
    logging.info("Checking if all specified files exist.")

    for file in files:
        if not os.path.isfile(file):
            logging.warning(f"File {file} does not exist.")
            return False

    logging.info("All specified files exist.")
    return True

def cansnp_file_writer(cansnp_data: List[List[Union[str, int]]], output_filepath: str) -> None:
    """
    Write Canonical SNPs (canSNPs) data to an output file.

    Args:
    cansnp_data: List of canSNP data, each represented as a list of strings or integers.
    output_filepath: Path to the output file where the canSNP data will be saved.

    Returns:
    None
    """
    try:
        logging.info(f"Starting to write canSNP data to {output_filepath}...")

        with open(output_filepath, 'w') as out_file:
            for i, cansnp_entry in enumerate(cansnp_data):
                formatted_entry = " ".join(map(str, cansnp_entry))
                out_file.write(f"{formatted_entry}\n")

                # Log progress every 10% if the data size is large enough
                if len(cansnp_data) > 10000:
                    current_percent = int(((i + 1) / len(cansnp_data)) * 100)
                    if current_percent % 10 == 0:
                        logging.info(f"Writing progress: {current_percent}%")

        logging.info(f"Successfully wrote {len(cansnp_data)} canSNP entries to {output_filepath}.")

    except Exception as e:
        logging.error(f"Failed to write canSNP data to {output_filepath}. Error: {str(e)}")
