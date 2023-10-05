import fasta_parser, yaml, re
import hc_parser
import sys, os.path
def read(file):
	with open(file,'r') as yam:
		return next(yaml.load_all(yam, Loader=yaml.FullLoader))

def clean_name(name):
	return re.sub('[^A-Za-z0-9]+', '', name)

def snp_extractor(fasta):
	alignment = list(fasta.values())	
	reference = alignment[0]
	alignment = alignment[1:]
	snps = []

	nucleotides = ['A','T','C','G']
	print("Extracting snps...")
	for index, nucleotide in enumerate(reference):
		snp = False
		if nucleotide not in nucleotides:
			continue
		variants = {nucleotide:1}
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
			snps.append([index,variants])
	print("Extraction of %i SNPs complete" %len(snps))	
	return snps

def canSNP_extractor(alignment, groups, snps):
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
		cansnps.append([str(variant), cansnp]+list(frequencies.keys()))
					
	print("Extraction of canSNPs complete...")
	return cansnps
def file_writer(data, output):
	with open(output, 'w') as out_file:
		for line in data:
			info = " ".join(line)
			out_file.write(info+"\n")
csvfile = sys.argv[1]
hcfile = sys.argv[2]
fastafile = sys.argv[3]
output = sys.argv[4]
if not os.path.isfile(csvfile) or not os.path.isfile(hcfile) or not os.path.isfile(fastafile) or os.path.isfile(output):
	print("One of the input files was not an existing file, exiting...")
	sys.exit(1) 
parser = hc_parser.hc_parser(csvfile)
yaml = read(hcfile)
hc = parser.reader(yaml)
fasta = fasta_parser.parse(fastafile, True, True)
snps = snp_extractor(fasta)
cansnps = canSNP_extractor(fasta,hc,snps)
file_writer(cansnps,output)
