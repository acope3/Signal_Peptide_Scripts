#!/bin/bash

for ((i=1;i<=500;i+=1))
do
	#perl ~/tai/misc/codonM Simulated/Simulated_nosp/nosp_sim_$i.fasta Simulated/Simulated_nosp/nosp_sim_$i.m
	perl ~/tai/misc/codonM Simulated/Simulated_nosp_aa_norm/nosp_sim_$i.fasta Simulated/Simulated_nosp_aa_norm/nosp_sim_$i.m
	#perl ~/tai/misc/codonM Simulated/Simulated_psp/psp_sim_$i.fasta Simulated/Simulated_psp/psp_sim_$i.m
	perl ~/tai/misc/codonM Simulated/Simulated_psp_aa_norm/psp_sim_$i.fasta Simulated/Simulated_psp_aa_norm/psp_sim_$i.m
	#perl ~/tai/misc/codonM Simulated/Simulated_sp/sp_sim_$i.fasta Simulated/Simulated_sp/sp_sim_$i.m
done
