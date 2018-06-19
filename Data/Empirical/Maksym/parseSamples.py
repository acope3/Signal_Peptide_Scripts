import re

gene_expression={}
with open("GSE45443_Transcripts_Samples_7_to_22.txt") as fin:
	for line in fin:
		line_spt=line.strip("\n").split("\t")
		if line_spt[5] != "-":
			gene_expression[line_spt[6]] = line_spt[8:]

with open("../../../Genomes/Ecoli_K12_ncbi/main_ht/mod_Ecoli_K12_main_ht_liberal.fasta",'r') as fin, open("samples_7_22.csv",'w') as out:
	out.write("Gene,Sample_7,Sample_8,Sample_9,Sample_10,Sample_11,Sample_12,Sample_13,Sample_14,Sample_15,Sample_16,Sample_17,Sample_18,Sample_19,Sample_20,Sample_21,Sample_22\n")
	#out.write("Gene,Sample_1,Sample_2,Sample_3,Sample_4,Sample_5,Sample_6\n")
	pat=re.compile("\[locus_tag=(b[0-9]{4})\]")
	for line in fin:
		if line[0] == ">":
			result = pat.search(line)
			out.write(line.split()[0][1:] + ",")
			if result == None:
				out.write("0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0\n")
			else:
				locus = result.group(1)
				exp = gene_expression.get(locus)
				if exp != None:
					out.write(",".join(exp))
					out.write("\n")
				else:
					out.write("0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0\n")