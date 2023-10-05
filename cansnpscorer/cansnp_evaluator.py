import sys, csv
import math, os.path
from operator import itemgetter
from codon_table import codon_table

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


#Method for checking distance to other snps
def neighbour_distance(snp, left, right):
	min_dist = 25
	left_dist = min(abs(left-snp),min_dist)
	right_dist = min(abs(right-snp),min_dist) 
	max_dist = min(left_dist, right_dist)
	return max_dist/float(min_dist)

#Calculate gc content
def gc_content(snp, reference, goal):
	distance = 50
	sequence = reference[snp-distance:snp+distance]
	GC = 0
	for i in sequence:
		if i == "C" or i == "G":
			GC += 1
	gc_content = GC/float(len(sequence))
	return 1-(abs(gc_content-goal)*(1/goal))

#Binary search algorithm
def search(snp, intervals, limit):
	found = False
	lowerbound = 0
	upperbound = len(intervals)-1
	while not found:
		if upperbound < lowerbound:
			return 1,	
		current = int(lowerbound + (upperbound-lowerbound)/2) 
		temp = intervals[current]

		left = max(temp[0]-limit,0)
		right = temp[1]+limit
		if snp >= left and snp <= right:
			if snp < temp[0]:
				return (temp[0]-snp)/float(limit),
			elif snp > temp [1]:
				return (snp-temp[1])/float(limit),
			else:
				return 0,temp
		elif snp < left:
			upperbound = current - 1
		elif snp > right:
			lowerbound = current + 1

#Method for checking type of substitution mutation	
def gene_checker(snp, reference, gene_interval, correction):
	variants = []
	index = snp[0] - gene_interval[0] - correction
	gene = reference[0][gene_interval[0]+correction:gene_interval[1]+correction+1]
	if gene_interval[2] == "-":
		index = (gene_interval[1] - gene_interval[0])- index
		gene = gene[::-1]

	j = index%3
	for i in snp[2:]:
		left = gene[index-j:index]
		mid = i
		right = gene[index+1:index+1+2-j]

		variants.append(left+mid+right)
	return codon_comparer(variants, gene_interval[2])
#Compare two codons
def codon_comparer(variants,strand):
		codon1 = variants[0]
		codon2 = variants[1]

		if strand == "-":
			codon1 = codon_complement(codon1)
			codon2 = codon_complement(codon2)
		N1 = codon_table[codon1]			
		N2 = codon_table[codon2]
		if N1 == N2:
			return 1
		elif N1 == "*" or N2 == "*":
			return 0
		else:
			return 0.5
#Translates codon in to their complement
def codon_complement(codon):
	translation = {"A":"T", "T":"A", "C":"G", "G":"C"}
	translated = ""
	for i in codon:
		translated = translated + translation[i]
	return translated
#Writes a csv file with snps
def file_writer(data, output):
	with open(output, 'w') as out_file:
		
		writer = csv.writer(out_file, delimiter = ',')
		writer.writerow(["SNP", "Group", "Score", "Neighbour-snps", "GC-content","tandem", "Gene", "Mutation"]) 
		for line in data:
			writer.writerow(line)

if len(sys.argv) != 6:
	print("Please enter a correct number of arguments...exiting")
	sys.exit(1) 


snps = snp_file_parser(sys.argv[1])
reference = reference_parser(sys.argv[2])
ref_gb = gff_parser(sys.argv[3])
trf = trf_parser(sys.argv[4])
output = sys.argv[5]

if os.path.isfile(output):
	print("Outputfile %s already exists... exiting")
	sys.exit(1)


modifier = 100

tandem_limit = 50
gene_limit = 0
gc_goal = 0.5

correction = 0
chromosome = 1

chromosome_length = reference[1]
scores = []
for i in range(0,len(snps)):
	if snps[i][1] == "None":
		continue	

	if snps[i][0] - correction > chromosome_length:
		correction = reference[chromosome]
		chromosome += 1
	if i > 0 and i < len(snps) - 1:
		neighbour = neighbour_distance(snps[i][0], snps[i-1][0], snps[i+1][0])

	elif i == 0:
		neighbour = neighbour_distance(snps[i][0], 0, snps[i+1])
	
	elif i == len(snps):
		neighbour = neighbour_distance(snps[i][0], snps[i-1], len(reference[0]))
	gc = gc_content(snps[i][0], reference[0], gc_goal)
	snp_pos = snps[i][0] - correction
	tandem = search(snp_pos, trf[chromosome], tandem_limit)
	gene = search(snp_pos, ref_gb[chromosome], gene_limit)
	mutation = 0

	if gene[0] == 0:
		mutation = gene_checker(snps[i],reference, gene[1], correction)
	score = [neighbour*1, gc*1 ,tandem[0]*1, (1 - gene[0])*2, mutation*1]
	score = [int(x*modifier) for x in score]
	total =  sum(score)
	
	scores.append([snps[i][0],snps[i][1], total] + score)

scores = sorted(scores, key=itemgetter(2), reverse = True)
file_writer(scores,output)
