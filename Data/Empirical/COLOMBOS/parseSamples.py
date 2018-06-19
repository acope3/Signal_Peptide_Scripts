import re

gene_expression={}
with open("GSE48829_counts.txt") as fin:
	for line in fin:
		line_spt=line.strip("\n").split("\t")
		gene_expression[line_spt[0]] = line_spt[1:]

with open("../../../Genomes/Ecoli_K12_ncbi/main_ht/mod_Ecoli_K12_main_ht_liberal.fasta",'r') as fin, open("samples.csv",'w') as out:
	#out.write("Gene,Sample_7,Sample_8,Sample_9,Sample_10,Sample_11,Sample_12,Sample_13,Sample_14,Sample_15,Sample_16,Sample_17,Sample_18,Sample_19,Sample_20,Sample_21,Sample_22\n")
	out.write("Gene,23A1,23B1,23C3,39A1,39B1,39C1,WTA1,WTB1,WTC3\n")
	pat=re.compile("\[locus_tag=(b[0-9]{4})\]")
	total = 0
	for line in fin:
		if line[0] == ">":
			result = pat.search(line)
			out.write(line.split()[0][1:] + ",")
			if result == None:
				out.write("0,0,0,0,0,0,0,0,0\n")
			else:
				locus = result.group(1)
				exp = gene_expression.get(locus)
				if exp != None:
					out.write(",".join(exp))
					out.write("\n")
				else:
					out.write("0,0,0,0,0,0,0,0,0\n")
					total +=1