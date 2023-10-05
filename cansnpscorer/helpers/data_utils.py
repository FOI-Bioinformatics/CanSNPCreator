from typing import List, Dict, Union, NoReturn

codon_table = {	"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
		"TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",	
		"TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",	
		"TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",	
		"CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
		"CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
		"CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
		"CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
		"ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
		"ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
		"AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
		"AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",	
		"GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",	
		"GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",	
		"GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",	
		"GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G" }

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

def clean_name(name: str) -> str:
	return re.sub('[^A-Za-z0-9]+', '', name)


def snp_extractor(fasta: Dict) -> List[Union[int, Dict]]:
	alignment = list(fasta.values())
	reference = alignment[0]
	alignment = alignment[1:]
	snps = []

	nucleotides = ['A', 'T', 'C', 'G']
	print("Extracting snps...")
	for index, nucleotide in enumerate(reference):
		snp = False
		if nucleotide not in nucleotides:
			continue
		variants = {nucleotide: 1}
		for genome in alignment:
			pos = genome[index]
			if pos not in variants:
				variants[pos] = 1
			else:
				variants[pos] += 1
			if pos not in nucleotides:
				snp = False
				break
			if pos != nucleotide:
				snp = True

		if snp and len(variants) == 2:
			snps.append([index, variants])
	print("Extraction of %i SNPs complete" % len(snps))
	return snps


def canSNP_extractor(alignment: Dict, groups: Dict, snps: List) -> List[List]:
	print("Extracting canSNPs...")
	cansnps = []
	for snp in snps:
		variant = snp[0]
		frequencies = snp[1]
		cansnp = "None"

		for group in groups:
			nucleotide = ""
			frequency = 1
			for index, isolate in enumerate(groups[group]):
				current = alignment[clean_name(isolate)][variant]
				if index == 0:
					nucleotide = current
				else:
					if nucleotide == current:
						frequency += 1
					else:
						break
			else:
				if frequency == frequencies[nucleotide]:
					cansnp = group
		cansnps.append([str(variant), cansnp] + list(frequencies.keys()))

	print("Extraction of canSNPs complete...")
	return cansnps