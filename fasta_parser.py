import re
import os.path
#Author Emil Samuelsson (2017-03-19)
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


