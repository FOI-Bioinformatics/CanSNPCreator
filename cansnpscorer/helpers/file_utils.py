
import re
import os.path
import sys
import csv
import yaml
from typing import List, Dict, Union, NoReturn


def clean_name(name):
	return re.sub('[^A-Za-z0-9]+', '', name)

#Two optional arguments names and dictonary. Supply names to get a list of list ([name,genome]).
#Use dictonary = True if you want a dictonary instead of list (requires names = True).
def fasta_parser(fasta, names, dictonary):
	genomes = []
	if dictonary:
		genomes = {}
	with open(fasta, 'r') as fast:
		open_fasta = fast.read().split(">")
		open_fasta = open_fasta[1:]
		for line in open_fasta:
			index = line.find("\n")
			name = clean_name(line[0 : index])
			genome = line[index + 1:].replace("\n","")
			if names:
				if dictonary:
					genomes[name] = genome
				else:
					genomes.append([name,genome])
			else:
				genomes.append(genome)
	return genomes


def parse(fasta, names = False, dictonary = False):
	if not os.path.isfile(fasta):
		print("File input does not exist. Enter a correct filepath")
		return -1
	if isinstance(names,bool) and isinstance(dictonary, bool):
		return fasta_parser(fasta, names, dictonary)
	else:
		print("Second and third arguments must be booleans")
		return -1


def clean_name(name):
    return re.sub('[^A-Za-z0-9]+', '', name)


class hc_parser():
    def __init__(self, csvfile):
        self.dictonary = {}
        self.group_label = "Group"
        self.group_index = 0
        self.strain_label = "Strain"
        self.strain_index = 0
        self.isolates = {}
        self.group_handler(csvfile)

    def group_handler(self, csv_file):
        with open(csv_file, 'r') as csvfile:
            reader = csv.reader(csvfile)
            for index, row in enumerate(reader):
                if index == 0:
                    if self.group_label in row and self.strain_label in row:
                        self.group_index = row.index(self.group_label)
                        self.strain_index = row.index(self.strain_label)
                    else:
                        print("Group or Strain column not found in csv file. Exiting...")
                        sys.exit(1)
                else:
                    if row[self.group_index] in self.isolates:
                        self.isolates[row[self.group_index]].append(clean_name(row[self.strain_index]))
                    else:
                        self.isolates[row[self.group_index]] = [clean_name(row[self.strain_index])]

    def dictonary_reader(self, yaml):
        isolates = []
        name = list(yaml)[0]
        items = list(yaml.values())[0]
        for group in items:
            if isinstance(group, dict):
                _isolates = self.dictonary_reader(group)
                isolates = isolates + _isolates
            else:
                group_isolate = self.isolates[group]
                isolates = isolates + group_isolate
                self.dictonary[group] = group_isolate
        self.dictonary[name] = isolates
        return isolates

    def reader(self, lista):
        for item in lista:
            if isinstance(item, dict):
                self.dictonary_reader(item)
            else:
                self.dictonary[item] = self.isolates[item]
        return self.dictonary


def snp_file_parser(snpfile):
	snps = []
	with open(snpfile,'r') as snp_file:
		allsnps = snp_file.read().split("\n")
		del allsnps[-1] #last element is empty due to split
	for snp in allsnps:
		if "#" in snp:
			continue
		else:
			try:
				variants = snp.split(" ")
				variants[0] = int(variants[0])
				snps.append(variants)
			except ValueError:
				print("ValueError, could not convert snp %s to an integer. Ignoring snp." %snp)
	return snps
def reference_parser(reference):
	ref = ""
	chromosomes = 0
	chro_lengths = []
	with open(reference,'r') as fil:
		temp = fil.read().split("\n")
		for line in temp:
			if ">" not in line:
				ref = ref + line.strip("\n")
			else:
				if chromosomes != 0:
					chro_lengths.append(len(ref))
				chromosomes += 1
	return_list = [ref] + chro_lengths
	return return_list

def gff_parser(gff_file):
	genes = {}
	with open(gff_file,'r') as gff:
		current_region = 0
		for line in gff:
			if "#" in line:
				continue
			gene = line.split("\t")
			if "region" in gene:
				current_region += 1
				genes[current_region] = []
			elif "gene" in gene:
				index = gene.index("gene")
				genes[current_region].append([int(gene[index+1])-1, int(gene[index+2])-1,gene[index+4]])
	return genes

def trf_parser(input):
	intervals = {}
	with open(input, 'r') as trf:
		chromosome = 0
		for line in trf:
			if "@" in line:
				chromosome += 1
				intervals[chromosome] = []
				continue

			temp_line = line.strip("\n")
			temp_line = temp_line.split(" ")
			try:
				first = int(temp_line[0])
				second = int(temp_line[1])
				intervals[chromosome].append([first,second])
			except ValueError:
				continue
	for i in intervals:
		intervals[i] = sorted(intervals[i], key=itemgetter(0))

	return intervals


# Writes a csv file with snps
def file_writer(data, output):
	with open(output, 'w') as out_file:
		writer = csv.writer(out_file, delimiter=',')
		writer.writerow(["SNP", "Group", "Score", "Neighbour-snps", "GC-content", "tandem", "Gene", "Mutation"])
		for line in data:
			writer.writerow(line)

def read(file: str) -> Dict:
	with open(file,'r') as yam:
		return next(yaml.load_all(yam, Loader=yaml.FullLoader))

def check_files_exist(files: List[str]) -> bool:
    return all(os.path.isfile(file) for file in files)

def cansnp_file_writer(data: List, output: str) -> None:
	with open(output, 'w') as out_file:
		for line in data:
			info = " ".join(line)
			out_file.write(info+"\n")