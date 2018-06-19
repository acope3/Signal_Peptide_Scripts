import pandas as pd
import numpy as np


transitions = ["GAA","AAA","TTA","GCC","TGC","GAC","TTC","GGC","CAC","CAA","ATC","CTG","AAC","CCC","CGC","TCC","ACC","GTC","TAC","AGC"]

def calc_mu(delta_M,neg=False):
	x = 0.5
	if neg:
		x *= -1
	mu = 1.0*10**(-9) * np.exp(x*delta_M)
	return mu


def calc_C0(average,total):
	c0 = (average*61.0 + total)/61.0
	return c0


def adj_dEta(df,aa):
	offset = df.loc[df["AA"] == aa,"Posterior"].min()
	if offset < 0.0:
		df.loc[df["AA"] == aa, "Posterior"] = df.loc[df["AA"] == aa,"Posterior"] - offset


def main():
	# df = pd.read_table("selection_mod_Ecoli_K12_MG1655_ncbi_main_liberal.csv",sep=",",header=0,usecols=[0,1,2])
	# aa_codon = pd.read_table("S.cerevisiae.tRNA.tsv",sep="\t",header=None,names=["AA","Codon"],usecols=[0,1])
	# df2 = pd.merge(left=aa_codon,right=df,how="left",on=["AA","Codon"])
	# df2.fillna(0,inplace=True)
	# aa = df2["AA"].unique()
	# for a in aa:
	# 	adj_dEta(df2,a)	
	
	
	# c0 = calc_C0(10.0,df2["Posterior"].sum())
	# #print c0
	# df3 = c0-df2.loc[:,"Posterior"] 
	# #print df3.mean()
	# df2.loc[:,"Posterior"] = df3
	# df2.to_csv("E.coli.tRNA.tsv",sep="\t",header=False,index=False)

	df = pd.read_table("mutation_mod_Ecoli_K12_MG1655_ncbi_main_liberal.csv",sep=",",header=0,usecols=[0,1,2])
	#aa_codon = pd.read_table("S.cerevisiae.tRNA.tsv",sep="\t",header=None,names=["AA","Codon"],usecols=[0,1])
	#df2 = pd.merge(left=aa_codon,right=df,how="left",on=["AA","Codon"])
	#df2.fillna(0,inplace=True)
	
	mu = calc_mu(df[["Posterior"]])
	mu2 = calc_mu(df[["Posterior"]],neg=True)
	df2 = pd.concat([df,mu],axis=1)
	df3 = pd.concat([df2,mu2],axis=1)
	transition_rate = 	0.0
	transversion_rate = 0.0
	total_trans = 0
	total_transv = 0
	for i in range(len(df2[["Codon"]])):
		if df3.iloc[i,1] in transitions:
			transition_rate = df3.iloc[i,3] + df3.iloc[i,4] + transition_rate
			total_trans +=2
		elif df3.iloc[i,1] not in ["CTT","CTC","CTA","AGA","AGA"]:
			transversion_rate = df3.iloc[i,3] + df3.iloc[i,4] + transversion_rate
			total_transv+=2
	print (transition_rate/total_trans)/(transversion_rate/total_transv)
	print total_trans,total_transv

main()