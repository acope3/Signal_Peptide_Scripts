## Author: Alex Cope
## Run in Python 2.7
## A script for performing acceptance-rejection sampling to generate the pseudo-secretory proteins
## discussed in the paper.
## File paths are hardcoded.
## Does allow the user some flexibility in how they perform the acceptance-rejection sampling.
## User may use using lognormal distribution (for \phi) or a normal distribution (for log(\phi)).
## Also have the choice of skewnorm distribution in case the log(\phi$) values are heavily skewed.
## Email Alex Cope at alexandercope0@gmail.com or acope3@vols.utk.edu

import numpy as np
import pandas as pd
import math
import random
from scipy.stats import skewnorm,norm,lognorm
from scipy.stats.mstats import gmean
import matplotlib.pyplot as plt
import sys 
seq ={}
genes = []
selected = []
phi = []
with open("../Data/Genomes/nosp_main.fasta",'r') as fin:
	gene = None
	for line in fin:
		if line[0] == ">":
			gene = line
			genes.append(gene)
		else:
			seq[gene] = line

if sys.argv[1] != "lognorm":
	exp_mc = pd.read_table("../Data/Empirical/Emp_Phi_Final/Mature_peptides/mp_main_phi.csv",sep=",",header=0)
	avg_mc = np.log(exp_mc[["Phi"]]).mean()


	exp_nosp = pd.read_table("../Data/Empirical/Emp_Phi_Final/nosp_main_phi.csv",sep=",",header=0)
	avg_nosp = np.log(exp_nosp[["Phi"]]).mean()
	

	mp = np.log(np.array(exp_mc[["Phi"]])[:,0])
	nosp = np.log(np.array(exp_nosp[["Phi"]])[:,0])
	print "Average log(Phi), secretory: ",avg_mc.values[0]
	print "Median log(Phi), secretory: ",np.median(mp)
	print "Variance log(Phi), secretory ", np.var(mp)
	print ""
	print "Average log(Phi), nonsecretory: ",avg_nosp.values[0]
	print "Median log(Phi), nonsecretory: ",np.median(nosp)
	print "Variance log(Phi), nonsecretory ", np.var(nosp)

else:
	exp_mc = pd.read_table("../Data/Empirical/Emp_Phi_Final/Mature_peptides/mp_main_phi.csv",sep=",",header=0)
	avg_mc = gmean(exp_mc[["Phi"]]) ##look at geometric mean because \phi follows a lognormal distribution, arithmetic mean will be biased

	exp_nosp = pd.read_table("../Data/Empirical/Emp_Phi_Final/nosp_main_phi.csv",sep=",",header=0)
	avg_nosp = gmean(exp_nosp[["Phi"]])

	mp = np.array(exp_mc[["Phi"]])[:,0]
	nosp = np.array(exp_nosp[["Phi"]])[:,0]
	print "Average Phi, secretory: ",avg_mc[0]
	print "Median Phi, secretory: ",np.median(mp)
	print "Variance Phi, secretory ", np.var(mp)
	print ""
	print "Average Phi, nonsecretory: ",avg_nosp[0]
	print "Median Phi, nonsecretory: ",np.median(nosp)
	print "Variance Phi, nonsecretory ", np.var(nosp)
total = 0.0



if sys.argv[1] == "skewnorm":
	mp_shape, mp_loc, mp_s,  = skewnorm.fit(mp)
	nosp_shape,nosp_loc, nosp_s = skewnorm.fit(nosp)
	dist_mp = skewnorm(mp_shape,mp_loc,mp_s)
	dist_nosp = skewnorm(nosp_shape,nosp_loc,nosp_s)
elif sys.argv[1] == "norm":
	mp_loc, mp_s = norm.fit(mp)
	nosp_loc, nosp_s = norm.fit(nosp)
	dist_mp = norm(mp_loc,mp_s)
	dist_nosp = norm(nosp_loc,nosp_s)
elif sys.argv[1] == "lognorm":
	mp_shape, mp_loc, mp_s,  = lognorm.fit(mp,floc=0)
	nosp_shape,nosp_loc, nosp_s = lognorm.fit(nosp,floc=0)
	dist_mp = lognorm(mp_shape,mp_loc,mp_s)
	dist_nosp = lognorm(nosp_shape,nosp_loc,nosp_s)


index = 0
for j in exp_nosp.iloc[:,1]:
	if sys.argv[1] != "lognorm":
		i= math.log(j)
	else:
		i = j
	pdf_mp = dist_mp.pdf(i)
	pdf_nosp = dist_nosp.pdf(i)
	c = float(sys.argv[2]) #specify c in the command line call, set so the assert statement below is true
	assert pdf_mp < (c*pdf_nosp)
	rand = random.uniform(0.0,1.0)
	if rand <= pdf_mp/(c*pdf_nosp):
		if len(seq.get(genes[index]).strip()[69:]) >= 45: #only include genes that have at least 15 codons after the first 23 codons. This choice is arbitrary
		#if len(seq.get(genes[index]).strip()) >= 300:
			selected.append(genes[index])
			phi.append(i)
			total += i
	index+=1
if sys.argv[1] == "skewnorm":
	shape,mu,s = skewnorm.fit(phi)
	print " "
	print "Shape Parameter: ", mp_shape,shape
	print "Location Parameter: ",mp_loc,mu
	print "Scale Parameter: ",mp_s,s
	print " "
	rv = skewnorm(shape,mu,s)
	print "Average",total/len(phi)
	print "Median", np.median(phi)
	print "Variance",np.var(phi)
	print " "
	print "Number",len(phi)
elif sys.argv[1] == "lognorm":
	shape,mu,s = lognorm.fit(phi,floc=0)
	print " "
	print "Shape Parameter: ", mp_shape,shape
	print "Location Parameter: ",mp_loc,mu
	print "Scale Parameter: ",mp_s,s
	print " "
	rv = lognorm(shape,mu,s)
	print "Average",np.mean(np.log(phi))
	print "Expected",np.log(exp_mc[["Phi"]]).mean()
	print "Median", np.median(phi)
	print "Variance",np.var(phi)
	print " "
	print "Number",len(phi)
else:
	mu,s = norm.fit(phi)
	print " "
	print "Location Parameter: ",mp_loc,mu
	print "Scale Parameter: ",mp_s,s
	print " "
	rv = norm(mu,s)
	print "Average",total/len(phi)
	print "Median", np.median(phi)
	print "Variance",np.var(phi)
	print " "
	print "Number",len(phi)

with open("../Data/Genomes/pseudo_sp_genes_w_emp_lognorm_final.fasta",'w') as full, open("../Data/Genomes/Pseudo_signal/pseudo_sp_w_emp_lognorm_final.fasta",'w') as sp, open("../Data/Genomes/Pseudo_mature/pseudo_mp_w_emp_lognorm_final.fasta",'w') as mp:
	for gene in selected:
		full.write(''.join([gene,seq.get(gene)]))
		sp.write(''.join([gene,seq.get(gene)[:69],'\n']))
		mp.write(''.join([gene,seq.get(gene)[69:]]))
assert len(selected) == len(phi)
with open("../Data/Empirical/Emp_Phi_Final/Pseudo_signal/pseudo_sp_w_emp_lognorm_final_phi.csv",'w') as out_1,open("../Data/Empirical/Emp_Phi_Final/Pseudo_mature/pseudo_mp_w_emp_lognorm_final_phi.csv",'w') as out_2:
	out_1.write("Gene_id,Phi\n")
	out_2.write("Gene_id,Phi\n")
	for i,j in zip(selected,phi):
		if sys.argv[1] != "lognorm":
			out_1.write(i.strip().split()[0][1:]+","+str(math.exp(j))+"\n")
			out_2.write(i.strip().split()[0][1:]+","+str(math.exp(j))+"\n")
		else:
			out_1.write(i.strip().split()[0][1:]+","+str(j)+"\n")
			out_2.write(i.strip().split()[0][1:]+","+str(j)+"\n")



