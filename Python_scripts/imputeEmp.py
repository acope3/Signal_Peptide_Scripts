## Author Alex Cope
## Run in Python 2.7
## Quick script for imputing missing values for fitting ROC-SEMPPR using empirical \phi values taken from Li et al 2014
## Values are imputed with the mean value of the empirical \phi values of those genes which have a ROC \phi value within 0.1
## of the ROC \phi value.
## This is based on the fact the ROC \phi values correlate well with the empirical \phi values 

import pandas as pd 
import numpy as np


ecoli_roc = pd.read_table("../Data/Empirical/ROC_Phi/Ecoli_K12_MG1655_main_phi.csv",sep=",",index_col=0,header=0)
li_emp = pd.read_table("../Data/Empirical/Li/li_mrna_level_main_ht_liberal.csv",sep=",",index_col=0,header=0).iloc[0:3358,:] ## only care about genes in the endogenous genome
eff = pd.read_table("../Data/Empirical/Li/li_trans_eff_main_ht_liberal.csv",sep=",",index_col=0,header=0).iloc[0:3358,:]
eff = eff.fillna(0)
li_emp = li_emp.assign(ppr=pd.Series(li_emp.iloc[:,0] * eff.iloc[:,0]).values)
ecoli_roc = ecoli_roc.apply(np.log)

sd = ecoli_roc.iloc[:,0].std(axis=0)
missing = li_emp[li_emp.iloc[:,1] == 0]
for gene in missing.index.values:
	roc_phi = ecoli_roc.loc[gene,:]
	similar_phi = ecoli_roc[abs(ecoli_roc.iloc[:,0]-roc_phi.values) < 0.1]
	emp_phi = li_emp.loc[similar_phi.index.values,:]
	emp_phi = emp_phi.loc[emp_phi.index.difference(missing.index)]
	avg_emp = emp_phi.iloc[:,0].mean()
	li_emp.loc[gene,"ppr"] = avg_emp
li_emp.to_csv("../Data/Empirical/Li/imputed_trans_eff.csv",sep=",",quoting=False,columns=["ppr"])