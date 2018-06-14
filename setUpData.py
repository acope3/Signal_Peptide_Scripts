## Author: Alex Cope
## Run in Python 2.7
## Command: python setUpData.py input.txt
## Contains many of the functions I used for setting up the data for analysis by ROC-SEMPPR
## This includes generating the signal peptide and mature peptide fasta files, mapping phi values
## to correpsonding genes, etc.
## Also added functionality for setting up data for CodonW analysis for determining "endogenous" and 
## "exogenous" genes. 
## Some functions are messier than others.
## Reads arguments, such as paths to appropriate files, from an input file, which is specified in the command line 

import sys
import os
import subprocess
import glob
import math
import re
import natsort


iso = ["ATT","ATC","ATA"]
leu = ["CTT","CTC","CTA","CTG","TTA","TTG"]
val = ["GTT","GTC","GTA","GTG"]
phe = ["TTT","TTC"]
cys = ["TGT","TGC"]
ala = ["GCT","GCC","GCA","GCG"]
gly = ["GGT","GGC","GGA","GGG"]
pro = ["CCT","CCC","CCA","CCG"]
thr = ["ACT","ACC","ACA","ACG"]
ser = ["TCT","TCC","TCA","TCG","AGT","AGC"]
tyr = ["TAT","TAC"]
trp = ["TGG"]
gln = ["CAA","CAG"]
asn = ["AAT","AAC"]
his = ["CAT","CAC"]
glu = ["GAA","GAG"]
asp = ["GAT","GAC"]
lys = ["AAA","AAG"]
arg = ["CGT","CGC","CGA","CGG","AGA","AGG"]
met = ["ATG"]
stop = ["TAA","TAG","TGA"]


amino_acids = [iso,leu,val,phe,cys,ala,gly,pro,thr,ser,tyr,trp,gln,asn,his,glu,asp,lys,arg,met,stop]

def readIDs(fileName):
	ids = {}
	with open(fileName,'r') as fin:
		fin.readline()
		for line in fin:
			line_spt = line.strip().split(",")
			ids[line_spt[0].split("_")[1]] = int(line_spt[5])
	return ids

def readGenes(fileName,regxp):
	loci = []
	pat = re.compile(regxp)
	with open(fileName) as fin:
		for line in fin:
			if line[0] == ">":
				results = pat.search(line)
				if results != None:
					locus = results.group(1)
				else:
					gene = None
				loci.append(locus)
	return loci

def writeGenes(loci,genes,fileName="loci_gene_mapping.csv"):
	assert len(loci) == len(genes)
	with open(fileName,'w') as out:
		out.write("Gene_Name,Locus_ID\n")
		for i in range(len(loci)):
			out.write(",".join([genes[i],loci[i]]))
			out.write("\n")

def totalPerGroup(ca_results,genome,regxp):
	ids = readIDs(ca_results)
	loci = readGenes(genome,regxp)
	totals = [0,0,0]
	group_3 =[]
	for i in range(0,len(loci)):
		value = ids.get(loci[i])
		if value == None:
			continue
		if value == 3:
			group_3.append(loci[i])
		ids.pop(loci[i])
		totals[value-1]+=1
	print totals
	return group_3


def filterGenome(exclude_list,fileName,regxp):
	pat = re.compile(regxp) 
	with open(fileName) as fin, open(fileName.split(".")[0]+"_main.fasta",'w') as out_1, open(fileName.split(".")[0]+"_ht.fasta",'w') as out_2:
		write_main = True
		write_limbo = False
		for line in fin:
			if line[0] == ">":
				results = pat.search(line)
				gene = results.group(1)
				if gene in exclude_list:
					write_main=False
					out_2.write(line)
				else:
					write_main = True
					out_1.write(line)
			else:
				if write_main:
					out_1.write(line)
				else:
					out_2.write(line)
	return None


def getTranslExcept(seq_desc):
	result = re.search("\[transl_except=\(pos:(\d+)\.\.(\d+),aa:\w+\)\]",seq_desc)
	return (int(result.group(1)),int(result.group(2)))

def removeUntranslatedRegion(seq,untranslated):
	seq_trim = "".join([seq[:untranslated[0]-1],seq[untranslated[1]:]])
	assert len(seq_trim.strip()) % 3 == 0
	return seq_trim


def forCodonW(genome,output,regxp):
	pat= re.compile(regxp)
	locus_seq = {}
	with open(genome) as fin, open(output,'w') as out:
		for line in fin:
			if line[0] == ">":
				total_nt = 0 
				is_transl_except = False
				result = pat.search(line)
				gene = result.group(1)
				seq_info = ">{}\n".format(gene)
				out.write(seq_info)
				if "transl_except" in line:
					untranslated_region = getTranslExcept(line)
					is_transl_except = True
			else:
				if not is_transl_except:
					out.write(line)
				else:
					seq_trim = removeUntranslatedRegion(line,untranslated_region)
					out.write(seq_trim)

	return None


def splitSeq(genome,output):
	with open(genome) as fin, open(output,'w') as out:
		for line in fin:
			if line[0] == ">":
				out.write(line)
			else:
				seq = line.strip()
				length = len(seq)
				total = 0
				segments = length/60
				for i in range(segments):
					#print total
					out.write(seq[total:total+60])
					out.write("\n")
					total = total + 60
				if total < length:
					out.write(seq[total:])
					out.write("\n")
	return None



def createCodonDict():
	codon_counts= {}
	codon_freq = {}
	codon_list = []
	with open("codons.txt") as fin:
		for line in fin:
			line_spt=line.split("\t")
			codon_counts[line_spt[1].strip()] = 0
			codon_freq[line_spt[1].strip()] = 0.0
			codon_list.append(line_spt[1].strip())
	return codon_counts, codon_freq, codon_list

def caiWeights():
	codon_counts, codon_freq, codon_list= createCodonDict()
	with open("ribo_prot.fasta",'r') as fin:
		for line in fin:
			pat = r"[ACGT]{3}"
			codons = re.findall(pat,line.strip())
			for i in codons:
				current = codon_counts.get(i)
				if current!=None:
					codon_counts[i] = current + 1
	with open("cai_weights_all.csv",'w') as out:
		out.write("Codon,Weight\n")
		for aa in amino_acids:
			counts = []
			for codon in aa:
				counts.append(codon_counts.get(codon))
			m = max(counts)
			for i in range(len(counts)):
				codon = aa[i]
				weight = float(counts[i])/m
				if counts[i] == 0:
					weight = 0.001
				out.write(codon +","+str(weight)+"\n")
def reorder():
	codon_dict = {}
	codon_array = [0 for i in range(64)]
	codon_order = ["" for i in range(64)]
	with open("codon_order.txt",'r') as fin:
		for line in fin:
			line_spt=line.strip().split()
			codon_dict[line_spt[1]] = int(line_spt[0])
	with open("cai_weights_all.csv",'r') as fin, open("codon_freq_codonw.txt",'w') as out:
		fin.readline()
		for line in fin:
			line_spt = line.strip().split(",")
			codon = line_spt[0]
			index = codon_dict.get(codon) - 1
			codon_array[index] = line_spt[1]
			codon_order[index] = codon
		for codon,weight in zip(codon_order,codon_array):
			out.write(str(weight)+" ") 

def geneSplit(gene,segments=300):
	ret = []
	length = len(gene)

	if length < 600:
		print "Gene is too short. Returning Empty list instead"
		return ret

	num_segments = length/segments

	end_length = length % segments
	start = 0
	end = segments
	for i in range(num_segments):
		ret.append(gene[start:end])
		start = end
		end = end + segments
	if end_length >= 45: #only keep segment if has 15 or more codons 
		ret.append(gene[start:])

	return ret


## out must end in "/"
def writeSegmentFasta(genes,out):
	cp ={} 
	cp.update(genes) ## Makes a copy of the dictionary in genes. Mutability of python dictionaries means following algorithm will eliminate data in genes, but might want to save
	numFiles = max(map(lambda x: len(x),genes.values())) ## Find the maximum gene length, know number of files to write
	files = []
	for i in range(1,numFiles+1):
		files.append(''.join([out,"segment_",str(i),".fasta"]))
	for i in range(numFiles):
		print files[i]
		with open(files[i],'w') as fin:
			for key in cp.keys():
				try:
					segment = cp.get(key)[i]
					assert len(segment) >= 45
					fin.write(''.join([key,segment,"\n"]))
				except IndexError:
					cp.pop(key)
	return None


def createDictionary(genome):
	genes = {}
	pseudo = False
	with open(genome,'r') as fin:
		gene_id = None
		for line in fin:
			if line[0] == ">":
				gene_id = line
			else:
				seg = geneSplit(line.strip(),75)
				if not seg:
					print "Skipping gene", gene_id
					continue
				genes[gene_id] = seg
	return genes 



def modFasta(inPath, outPath,sep="\n",index=0,include_pseudo=False):
	first = True
	pseudo = False
	with open(inPath,'r') as fin, open(outPath,'w') as out:
		for line in fin:
			if line[0] == ">":
				if "[pseudo=true]" in line:
					pseudo = True
					if not include_pseudo:
						continue
				if first:
					out.write(line.split(sep)[index])
					out.write("\n")
					first = False
				else:
					out.write("\n")
					out.write(line.split(sep)[index])
					out.write("\n")
				pseudo = False
			else:
				if pseudo and not include_pseudo:
					continue
				else:
					out.write(line.strip())
	return None

def readInGenome(inPath):
	genes = []
	with open(inPath,'r') as fin:
		for line in fin:
			if line[0] == ">":
				line_spt = line.split()
				genes.append(line_spt[0][1:])
	return genes

def getPhi(genes, inPath, ids = False, sep = ',', log = True, header = None,col=0):
	expression = {}
	index = 0
	with open(inPath,'r') as exp:
		if header:
			exp.readline()
		for line in exp:
			value = float(line.strip().split(sep)[col])
			if log:
				expression[genes[index]] = math.pow(10,value)
			else:
				expression[genes[index]] = value
			index+=1
	return expression

def createPhiFile(inPath,outPath,expression):
	files = glob.glob(inPath+"/*.fasta")
	for fil in files:
		fil_spt = fil.split("/")
		name = fil_spt[len(fil_spt)-1][:-6]
		with open(fil,'r') as fin, open(outPath+"/"+name+"_phi.csv",'w') as out:
			out.write("Gene_id,Phi\n")
			for line in fin:
				if line[0] == ">":
					key = line.split()[0][1:]
					phi = expression.get(key)
					if phi != None:
						out.write(key+","+str(phi)+"\n")
	return None

def getGenesWSigPep(inPath,score=0.5):
	genes = {}
	with open(inPath,'r') as fin:
		fin.readline()
		fin.readline()
		fin.readline()
		for line in fin:
			line = line.split()
			if float(line[5]) > score:
				genes[line[0]] = (int(line[4]),float(line[5]))
	return genes


def createGenomes(sp_genes,genome,nosp,sp,mp,genes_w_sp,regxp):
	cut_off = -1
	pat = re.compile(regxp)
	pseudo = False
	with open(genome,'r') as fin, open(sp,'w') as sp_fasta, \
	 open(mp,'w') as mp_fasta, open(nosp,'w') as nosp_fasta, \
	 open(genes_w_sp,'w') as sp_genes_fasta:
		value = None
		for line in fin:
			if line[0] == ">":
				try:
					value = sp_genes.get(pat.search(line).group())
				except:
					pseudo = True
					continue
				if value != None:
					print pat.search(line).group()
					cut_off = value[0] * 3
					sp_fasta.write(line)
					mp_fasta.write(line)
					sp_genes_fasta.write(line)
				else:
					nosp_fasta.write(line)
					cut_off = -1
			elif not pseudo:
				if cut_off != -1:
					sp_fasta.write(line[:cut_off])
					sp_fasta.write("\n")
					mp_fasta.write(line[cut_off:])
					sp_genes_fasta.write(line)
				else:
					nosp_fasta.write(line)
			else:
				pseudo = False

def readInputFile(inPath):
	values = {}
	with open(inPath,'r') as fin:
		for line in fin:
			line_spt = line.split(":")
			if line_spt[0] == "Genome":
				values["Genome"] = line_spt[1].strip()
			elif line_spt[0] == "Modified Genome":
				values["Modified Genome"] = line_spt[1].strip()
			elif line_spt[0] == "CodonW Genome":
				values["CodonW Genome"] = line_spt[1].strip()
			elif line_spt[0] == "SignalP Output":
				values["SignalP Output"] = line_spt[1].strip()
			elif line_spt[0] == "Genome Folder":
				values["Genome Folder"] = line_spt[1].strip()
			elif line_spt[0] == "Expression Input":
				values["Expression Input"] = line_spt[1].strip()
			elif line_spt[0] == "Expression Output":
				values["Expression Output"] = line_spt[1].strip()
			elif line_spt[0] == "Mature Peptide Fasta":
				values["Mature Peptide Fasta"] = line_spt[1].strip()
			elif line_spt[0] == "Signal Peptide Fasta":
				values["Signal Peptide Fasta"] = line_spt[1].strip()
			elif line_spt[0] == "Non-signal peptide Fasta":
				values["Non-signal peptide Fasta"] = line_spt[1].strip()
			elif line_spt[0] == "Signal peptide genes Fasta":
				values["Signal peptide genes Fasta"] = line_spt[1].strip()
			elif line_spt[0] == "Protein ID RegXP":
				values["Protein ID RegXP"] = line_spt[1].strip()
			elif line_spt[0] == "Locus ID RegXP":
				values["Locus ID RegXP"] = line_spt[1].strip()
			elif line_spt[0] == "Filter":
				values["Filter"] = line_spt[1].strip()

	return values



def prepForCodonW():
	values = readInputFile(sys.argv[1])
	genome = values["Genome"]
	mod_genome = values["Modified Genome"]
	modFasta(genome,mod_genome,include_pseudo=True)

	codonw = values["CodonW Genome"]
	regxp = values["Locus ID RegXP"]
	ca_results = values["Filter"]
	exclude_list=totalPerGroup(ca_results,genome,regxp)
	filterGenome(exclude_list,genome,regxp)
	splitSeq("tmp.txt","tmp_2.txt")
	forCodonW("tmp_2.txt",codonw,regxp)
	os.remove("tmp.txt")
	os.remove("tmp_2.txt")
	return 0


def prepForROCAnalyses():
	values = readInputFile(sys.argv[1])
	genome = values["Genome"]
	# mod_genome = values["Modified Genome"]
	# modFasta(genome,mod_genome,include_pseudo=False)
	# sp_out = values["SignalP Output"]
	# sp_genes = getGenesWSigPep(sp_out)
	# mp_fasta = values["Mature Peptide Fasta"]
	# sp_fasta = values["Signal Peptide Fasta"]
	# nosp_fasta = values["Non-signal peptide Fasta"]
	# sp_genes_fasta = values["Signal peptide genes Fasta"]
	# regxp = values["Protein ID RegXP"]
	# createGenomes(sp_genes,genome,nosp_fasta,sp_fasta,mp_fasta,sp_genes_fasta,regxp)
	genes = readInGenome(genome)
	xp_in = values["Expression Input"]
	expression = getPhi(genes,xp_in,log=False,header=True,col=1)
	gen_folder = values["Genome Folder"]
	#if not os.path.exists(gen_folder):
	#	os.makedirs(gen_folder)
	#out = createDictionary(mod_genome,600)
	#writeSegmentFasta(out,gen_folder)
	xp_out = values["Expression Output"]
	if not os.path.exists(xp_out):
		os.makedirs(xp_out)
	createPhiFile(gen_folder,xp_out,expression)
	return 0

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "Incorrect number of inputs. User should specify path to input file for "
	if sys.argv[2] == "CodonW":
		prepForCodonW()
	elif sys.argv[2] == "ROC":
		prepForROCAnalyses()
	else:
		print "Second input should be either CodonW or ROC"