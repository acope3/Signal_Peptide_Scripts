## Author: Alex Cope
## Run with Python 2.7
## An ugly script written to map gene expression data.
## File paths are currently hardcoded
## Also checks for gene synonyms (taken from bioDBnet) when performing the mapping
## If you have questions, do not hesitate to email me at alexandercope0@gmail.com or acope3@vols.utk.edu

import re

def getEmpPhi(fileInput,sep=",",metric_index=1,low_conf=True):
	mapping = {}
	with open(fileInput,'r') as fin:
		fin.readline() # read header
		for line in fin:
			line_spt = line.strip().split(sep)
			if "[" in line_spt[metric_index]:
				if low_conf:
					value = "0"
				else:
					value = line_spt[metric_index][1:-1]
			else:
				value = line_spt[metric_index] 
			mapping[line_spt[0]] = value 
	return mapping

def geneSyn(fileName):
	synoyms = {}
	with open(fileName,'r') as fin:
		fin.readline()
		for line in fin:
			line_spt = line.strip().split("\t")
			synoyms[line_spt[0]] = line_spt[1].split("; ")
	return synoyms

def main():
	synoyms = geneSyn("../Data/Empirical/Ecoli_K12_ncbi/Li/gene_syn_mapping.txt")
	mRNA_mapping = getEmpPhi("../Data/Empirical/Ecoli_K12_ncbi/Li/prot_synth.tsv",sep="\t",metric_index=3,low_conf=False)
	trans_eff = getEmpPhi("../Data/Empirical/Ecoli_K12_ncbi/Li/mrna.tsv",sep="\t",metric_index=2,low_conf=False)
	with open("../Data/Genomes/Ecoli_K12_ncbi/main_ht/mod_Ecoli_K12_main_ht_liberal.fasta",'r') as fin, open("../Data/Empirical/Ecoli_K12_ncbi/Li/li_prot_synth_main_ht_liberal_mops_complete_no_meth_low_conf_included.csv",'w') as mRNA_out, \
		open("../Data/Empirical/Ecoli_K12_ncbi/Li/li_trans_eff_main_ht_liberal.csv",'w') as trans_out:
		mRNA_out.write("Gene,mRNA_RPKM\n")
		trans_out.write("Gene,Translation_Efficiency_AU\n")
		total = 0
		ids = mRNA_mapping.keys()
		for line in fin:
			if line[0] == ">":
				if "[pseudo=true]" not in line:
					line_spt = line.split()
					mRNA_out.write(line_spt[0][1:] + ",")
					trans_out.write(line_spt[0][1:] + ",")
					gene = line_spt[1][6:-1]
					if "ins" in gene:
						gene = gene[:4]
					rpkm = mRNA_mapping.get(gene) #[gene=xxxX]
					au = trans_eff.get(gene)
					if rpkm == None:
						syns = synoyms.get(gene)
						tmp = None
						for i in syns:
							if i in ids:
								tmp = mRNA_mapping.get(i)
								mRNA_out.write(tmp+"\n")
								break
						if tmp == None:
							mRNA_out.write("0"+"\n")
							print gene
							total+=1
					else:
					 	mRNA_out.write(rpkm+"\n")
					if au == None:
						syns = synoyms.get(gene)
						tmp = None
						for i in syns:
							if i in ids:
								tmp = trans_eff.get(i)
								trans_out.write(tmp+"\n")
								break
						if tmp == None:
							trans_out.write("0.0"+"\n")
							print gene
							total+=1
					elif au == "NA":
					  	trans_out.write("0.0"+"\n")
					else:
					  	trans_out.write(au+"\n")
		print total

main()