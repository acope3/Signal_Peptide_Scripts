## Author: Alex Cope
## Email: acope3@vols.utk.edu, alexandercope0@gmail.com
## A script hacked together for checking frequencies of amino acids in signal peptides and 5'-ends of nonsecretory proteins and to generate
## random peptide sequences with amino acid frequencies comparable to signal peptides. 
## Not written well because I did not originally intend to put it on github, but I decided to and did not want to clean it up to much
## for the sake of reproducibility.
## Expects two commandline inputs, the file with the signal peptide sequences and either the 5'-end sequences or the pseudo-signal peptides
##

import sys
import re
import random 
from scipy.stats import chisquare,pearsonr
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pylab import savefig

iso = ["ATT","ATC","ATA"]
leu = ["CTT","CTC","CTA","CTG","TTA","TTG"]
val = ["GTT","GTC","GTA","GTG"]
phe = ["TTT","TTC"]
cys = ["TGT","TGC"]
ala = ["GCT","GCC","GCA","GCG"]
gly = ["GGT","GGC","GGA","GGG"]
pro = ["CCT","CCC","CCA","CCG"]
thr = ["ACT","ACC","ACA","ACG"]
ser1 = ["TCT","TCC","TCA","TCG"]
ser2= ["AGT","AGC"]
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


aa_dict = {"A":ala,
		   "C":cys,
		   "D":asp,
		   "E":glu,
		   "F":phe,
		   "G":gly,
		   "H":his,
		   "I":iso,
		   "K":lys,
		   "L":leu,
		   "M":met,
		   "N":asn,
		   "P":pro,
		   "Q":gln,
		   "R":arg,
		   "S":ser1,
		   "T":tyr,
		   "V":val,
		   "W":trp,
		   "Y":thr,
		   "Z":ser2}

def getAA(codon):
	for a in aa_dict.keys():
		syn_codons = aa_dict.get(a)
		if codon in syn_codons:
			return a


def calculateAAFreq(fileName):
	codon_freq = {}
	aa_freq = {}
	total_aa = 0
	pat = r"[ACGT]{3}"
	with open(fileName,'r') as fin:
		for line in fin:
			if line[0] != ">":
				seq = re.findall(pat,line.strip())
				#not counting methionine start codon from frequency, ROC-SEMPPR automatically adds this
				for codon in seq[1:]:
					if codon == "TAA" or codon == "TAG" or codon == "TGA":
						continue
					try:
						value = codon_freq.get(codon) + 1
						codon_freq[codon] = value
					except:
						codon_freq[codon] = 1.0
					total_aa += 1
		for aa in aa_dict.keys():
			total = 0
			for codon in aa_dict.get(aa):
				try:
					total = total + codon_freq.get(codon)
				except:
					continue
			aa_freq[aa] = total/total_aa
	return aa_freq,total_aa

def calculateAAFreq2(fileName):
	codon_freq = {}
	aa_freq = {}
	total_aa = 0
	pat = r"[ACGT]{3}"
	with open(fileName,'r') as fin:
		for line in fin:
			if line[0] != ">":
				seq = re.findall(pat,line.strip())
				assert len(seq) == 23
				for codon in seq[1:]:
					if codon == "TAA" or codon == "TAG" or codon == "TGA":
						continue
					try:
						value = codon_freq.get(codon) + 1
						codon_freq[codon] = value
					except:
						codon_freq[codon] = 1.0
					total_aa += 1
		for aa in aa_dict.keys():
			total = 0
			for codon in aa_dict.get(aa):
				try:
					total = total + codon_freq.get(codon)
				except:
					continue
			aa_freq[aa] = total
	return aa_freq,total_aa

def probRange(aa_freq,value):
	lower = 0.0
	amino = aa_freq.keys()
	random.shuffle(amino)
	for a in amino:
		if lower <= value and value < (lower + aa_freq.get(a)):
			return a
		else:
			lower += aa_freq.get(a)

def checkAAFreq(freq_dict):
	for aa in freq_dict.keys():
		print aa, freq_dict.get(aa)

sp_freq,sp_total = calculateAAFreq(sys.argv[1])
nosp_freq,nosp_total = calculateAAFreq(sys.argv[2])

# AA = np.repeat(aa_dict.keys(), 2)
# print AA
# freq = []
# for a in aa_dict.keys():
# 	freq.append(sp_freq.get(a)*100)
# 	freq.append(nosp_freq.get(a)*100)


# df = pd.DataFrame(data={"Amino acid":AA,"Region":np.tile(np.array(["Signal Peptide","5'-end (Nonsecretory Proteins)"]),21),"Frequency(%)":freq})
# a4_dims = (7.5, 5)
# fig, ax = plt.subplots(figsize=a4_dims)
# sns.barplot(ax=ax,x="Amino acid",y="Frequency(%)",hue="Region",data=df,ci=None)
# fig = ax.get_figure()
# fig.savefig('../aafreq.pdf', bbox_inches='tight',dpi=300)
#plt.show()
# checkAAFreq(sp_freq)
#checkAAFreq(nosp_freq)


#pat = re.compile(r'[AGCT]{3}')
cur_freq = {}
for a in aa_dict.keys():
	cur_freq[a] = 0
cur_total = 0.0
tries = 0
with open(sys.argv[2]) as fin, open("pseudo_sp_aa_norm.fasta",'w') as out:
	gene = ''
	for line in fin:
		if line[0] == ">":
			gene = line
		else:
			seq =re.findall(r'[AGCT]{3}',line.strip())
			## start sequence with start codon
			new_seq = [seq[0]]
			total = 1
			while total < 23:
				if cur_total == 0:
					rand = random.random()
					aa = probRange(sp_freq,rand)
				else:
					check = False
					while not check:
						rand = random.random()
						aa = probRange(sp_freq,rand)
						current = cur_freq.get(aa)
						current = current/cur_total
						if current <= sp_freq.get(aa):
							check = True

				occur_seq = []
				for j in seq[1:]:
					tmp = getAA(j)
					if tmp == aa:
						occur_seq.append(j)
				if len(occur_seq) != 0:
					index = random.randint(0,len(occur_seq)-1)
					new_seq.append(occur_seq[index])
					cur_freq[aa] = cur_freq.get(aa) + 1
					cur_total += 1
					total += 1
				else:
					# tries += 1
					# codons = aa_dict.get(aa)
					# index = random.randint(0,len(codons)-1)
					# new_seq.append(codons[index])
					# cur_freq[aa] = cur_freq.get(aa) + 1
					# cur_total += 1
					total += 1
			if len(new_seq) >=10:
				out.write(gene)
				out.write("".join(new_seq))
				out.write("\n")
print tries
new_psp, new_total = calculateAAFreq("pseudo_sp_aa_norm.fasta")
checkAAFreq(new_psp)	
print " "
checkAAFreq(sp_freq)	

