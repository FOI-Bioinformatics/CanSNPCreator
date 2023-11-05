
import re
import os.path
import sys
import csv
import logging
from ruamel.yaml import YAML
from ruamel.yaml.representer import Representer
from typing import List, Dict, Union, NoReturn, Tuple, Optional
from operator import itemgetter
import json
from collections import defaultdict
from io import StringIO

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

def extract_yaml_data(data, data_set=None):
    """
    Recursively extract all keys and values from a nested structure (dict or list)
    and store them in a set. The set is then sorted and returned as a list.

    Args:
    - data (Union[dict, list]): The input data structure.
    - data_set (set, optional): Set to store the keys and values. Defaults to an empty set.

    Returns:
    - list: A sorted list containing all unique keys and values from the data.
    """
    if data_set is None:
        data_set = set()

    if isinstance(data, dict):
        for key, value in data.items():
            data_set.add(key)
            extract_yaml_data(value, data_set)
    elif isinstance(data, list):
        for item in data:
            extract_yaml_data(item, data_set)
    else:
        data_set.add(data)

    return sorted(data_set)


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

    def __init__(self, group_data: str, data_format: str = 'csv') -> None:
        """
        Initialize the HCParser object.

        Args:
            group_data (str): The path to the file containing group and strain information.
            data_format (str): The format of the group data ('csv' or 'repgenr_group').
        """
        logging.info(f"Initializing HCParser with data file: {group_data}")
        self.dictionary: Dict[str, List[str]] = {}
        self.group_label = "Group"
        self.strain_label = "Strain"
        self.isolates: Dict[str, List[str]] = {}
        if data_format == 'csv':
            self._group_handler(group_data)
        elif data_format == 'repgenr_group':
            self._repgenr_group_handler(group_data)
        else:
            raise ValueError(f"Unknown data format: {data_format}")

    def _repgenr_group_handler(self, repgenr_group_data_file: str) -> None:
        """
        Process the repgenr_group data file and populate the 'isolates' dictionary.

        Args:
            repgenr_group_data_file (str): The path to the repgenr_group data file.
        """
        logging.info(f"Processing repgenr_group data file: {repgenr_group_data_file}")
        with open(repgenr_group_data_file, 'r') as file:
            reader = csv.DictReader(file, delimiter='\t')
            for row in reader:
                representative = row['representative'].replace('.fasta', '')
                contained = row['contained'].replace('.fasta', '')
                is_self = row['is_self'] == 'True'
                if is_self:
                    # If the row indicates that the strain is representative of itself, create a new group
                    self.isolates[representative] = [contained]
                else:
                    # Otherwise, add the strain to the group of its representative
                    self.isolates.setdefault(representative, []).append(contained)


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
    header = ["Position", "Group", "Score", "NeighborDistanceScore", "GCContentScore", "TandemRepeatProximityScore", "GeneProximityScore", "MutationTypeScore"]

    # Log the number of CanSNPs entries
    logging.info(f"Writing {len(snp_data)} CanSNP entries to the file {output_file_path}.")

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
                        logging.info(f"Writing canSNP data progress: {current_percent}%")

        logging.info(f"Successfully wrote {len(cansnp_data)} canSNP entries to {output_filepath}.")

    except Exception as e:
        logging.error(f"Failed to write canSNP data to {output_filepath}. Error: {str(e)}")

def recognize_group_file_format(file_path: str) -> str:
    """
    Recognize the file format by inspecting the first few lines.

    Args:
        file_path (str): The path to the file to inspect.

    Returns:
        str: The recognized file format ('csv' or 'repgenr_group').
    """
    logging.basicConfig(level=logging.INFO)
    try:
        with open(file_path, 'r') as file:
            lines: List[str] = [file.readline().strip() for _ in range(5)]  # Read the first 5 lines
            header: str = lines[0]  # Typically the header should be the first line

            # Define headers for the formats
            csv_header: str = 'Strain,Group'
            repgenr_group_header: str = 'representative\tcontained\tis_self'

            # Check if any of the first five lines starts with the expected headers
            if any(line.startswith(csv_header) for line in lines):
                logging.info(f"The file at '{file_path}' is recognized as CSV format.")
                return 'csv'
            elif any(line.startswith(repgenr_group_header) for line in lines):
                logging.info(f"The file at '{file_path}' is recognized as repgenr_group format.")
                return 'repgenr_group'
            else:
                raise ValueError(f"Unknown file format for file: {file_path}")

    except FileNotFoundError:
        logging.error(f"The file at '{file_path}' was not found.")
        raise
    except Exception as e:
        logging.error(f"An error occurred while trying to read the file: {e}")
        raise

def recognize_file_format(file_path: str) -> Union[str, None]:
    logging.info(f"Attempting to recognize the format of file: {file_path}")

    with open(file_path, 'r') as f:
        first_lines = [next(f).strip() for _ in range(10) if not f.readline().strip() == ""]

    if not first_lines:
        logging.warning("The file appears to be empty or only contains whitespace.")
        return None

    # Check if it's a Newick format
    if any(";" in line for line in first_lines):
        logging.info(f"The file {file_path} appears to be in Newick format.")
        return "newick"

    # Check if it's a JSON format
    if any(line.startswith("{") or line.startswith("[") for line in first_lines):
        logging.info(f"The file {file_path} appears to be in JSON format.")
        return "json"

    # Check if it's a YAML format
    if any(":" in line for line in first_lines):
        logging.info(f"The file {file_path} appears to be in YAML format.")
        return "yaml"

    logging.warning(f"Unable to determine the format of the file {file_path}.")
    return "unknown"


def convert_data(data: str, input_format: str, output_format: str) -> str:
    logging.info(f"Starting conversion from {input_format} to {output_format}")

    def hier_to_child_parent_with_leaf(hier: Union[Dict, List], parent=None, cp_list=[]):
        logging.debug(f"Converting hierarchy to child-parent: {hier} with parent {parent}")
        if isinstance(hier, dict):
            for k, v in hier.items():
                cp_list.append((k, parent))
                hier_to_child_parent_with_leaf(v, k, cp_list)
        elif isinstance(hier, list):
            for i in hier:
                hier_to_child_parent_with_leaf(i, parent, cp_list)
        else:
            cp_list.append((hier, parent))
        logging.debug(f"Child-parent list so far: {cp_list}")
        return cp_list

    def child_parent_to_hier_with_leaf(cp_list: List[Tuple[str, str]], root_is_list: bool = False) -> Union[Dict, List]:
        logging.debug(f"Converting child-parent list to hierarchy: {cp_list}")
        parent_to_children = defaultdict(list)

        for child, parent in cp_list:
            parent_to_children[parent].append(child)

        def build_hier_from_parent_with_leaf(parent: Optional[str]) -> Union[Dict, List]:
            children = parent_to_children.get(parent, [])
            if not children:
                return None
            if len(children) != len(set(children)):
                return [build_hier_from_parent_with_leaf(child) for child in children]
            else:
                return {child: build_hier_from_parent_with_leaf(child) for child in children}

        hier = build_hier_from_parent_with_leaf(None)
        logging.debug(f"Hierarchy so far: {hier}")

        return hier

        def hier_to_child_parent(hier: Union[Dict, List], parent=None, cp_list=[]):
            if isinstance(hier, dict):
                for k, v in hier.items():
                    cp_list.append((k, parent))
                    hier_to_child_parent(v, k, cp_list)
            elif isinstance(hier, list):
                for i in hier:
                    hier_to_child_parent(i, parent, cp_list)
            else:
                cp_list.append((hier, parent))
            return cp_list

        def child_parent_to_hier(cp_list: List[Tuple[str, str]], root_is_list: bool = False) -> Union[Dict, List]:
            parent_to_children = defaultdict(list)

            for child, parent in cp_list:
                parent_to_children[parent].append(child)

            def build_hier_from_parent(parent: str) -> Union[Dict, List]:
                children = parent_to_children.get(parent, [])
                if not children:
                    return None
                if len(children) != len(set(children)):
                    return [build_hier_from_parent(child) for child in children]
                else:
                    return {child: build_hier_from_parent(child) for child in children}

            if root_is_list:
                return [build_hier_from_parent(None)]
            else:
                return build_hier_from_parent(None)


        def build_hier_from_parent_with_leaf(parent: Optional[str]) -> Union[Dict, List]:
            children = parent_to_children.get(parent, [])
            if not children:
                return None
            if len(children) != len(set(children)):
                return [build_hier_from_parent_with_leaf(child) for child in children]
            else:
                return {child: build_hier_from_parent_with_leaf(child) for child in children}

        if root_is_list:
            return [build_hier_from_parent_with_leaf(None)]
        else:
            return build_hier_from_parent_with_leaf(None)

    def remove_leaf_colons(yaml_str: str) -> str:
        lines = yaml_str.splitlines()
        processed_lines = []
        for line in lines:
            if line.strip().endswith(":"):
                processed_lines.append(line.rstrip(":"))
            else:
                processed_lines.append(line)
        return "\n".join(processed_lines)

    def newick_to_child_parent(newick_str: str, parent=None, cp_list=[]):
        stack = []
        child_name = ""
        for char in newick_str:
            if char == "(":
                stack.append((child_name, parent))
                parent = child_name
                child_name = ""
            elif char == "," or char == ")":
                stack.append((child_name, parent))
                cp_list.append((child_name, parent))
                child_name = ""
            elif char != ";":
                child_name += char
        return cp_list

    def build_newick(node, cp_list):
        children = [child for child, parent in cp_list if parent == node]
        if children:
            return f"({','.join([build_newick(child, cp_list) for child in children])}){node}"
        else:
            return node

    # Custom representer for dictionaries with null values
    def null_dict_representer(dumper, data):
        return dumper.represent_dict({k: v for k, v in data.items() if v is not None})

    def remove_empty_dicts(yaml_str):
        return yaml_str.replace(": {}", "")

    # Function to dump YAML with custom representer
    # Function to dump YAML with custom representer
    def dump_custom_yaml(data):
        yaml = YAML()
        yaml.default_flow_style = None  # Use block style for collections
        yaml.indent(sequence=4, offset=2)  # Indentation settings
        yaml.representer.add_representer(dict, null_dict_representer)  # Add custom representer to this instance

        stream = StringIO()
        yaml.dump(data, stream)
        yaml_str = stream.getvalue()
        yaml_str = remove_empty_dicts(yaml_str)

        return yaml_str

    logging.info("Reading input data...")
    if input_format == 'yaml':
        yaml = YAML()
        input_data = yaml.load(data)
    elif input_format == 'json':
        input_data = json.loads(data)
    elif input_format == 'newick':
        input_data = newick_to_child_parent(data)
    elif input_format == 'child_parent':
        input_data = [(line.split()[0], line.split()[1] if len(line.split()) > 1 else None) for line in
                      data.strip().split('\n') if len(line.split()) > 1]
    else:
        logging.error("Unsupported input format")
        raise ValueError("Unsupported input format")

    logging.info("Converting to child-parent list...")
    if input_format in ['yaml', 'json']:
        cp_list = hier_to_child_parent_with_leaf(input_data)
    else:
        cp_list = input_data

    logging.info("Converting to output format...")
    if output_format == 'yaml':
        output_data = child_parent_to_hier_with_leaf(cp_list, root_is_list=isinstance(input_data, list))
        return dump_custom_yaml(output_data)
    elif output_format == 'json':
        output_data = child_parent_to_hier_with_leaf(cp_list, root_is_list=isinstance(input_data, list))
        return json.dumps(output_data, indent=2)
    elif output_format == 'newick':
        root = [node for node, parent in cp_list if parent is None][0]
        return build_newick(root, cp_list) + ';'
    elif output_format == 'child_parent':
        return "\n".join([f"{child}\t{parent}" for child, parent in cp_list if parent is not None])
    else:
        logging.error("Unsupported output format")
        raise ValueError("Unsupported output format")

def save_to_file(data: str, filename: str, file_format: str):
    """
    Save converted data to a file.

    Parameters:
        - data (str): The converted data as a string
        - filename (str): The name of the output file
        - file_format (str): The format of the data ("yaml", "json", "newick", "child_parent")

    """
    with open(f"{filename}.{file_format}", "w") as f:
        # Adding a header for child_parent format
        if file_format == "child_parent":
            f.write("child\tparent\n")
        f.write(data)
