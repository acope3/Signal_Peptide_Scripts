## Author: Alex Cope
## Script used for simulating codon usage in signal peptides, 5'-ends, and pseudo-signal peptides
## Need to specify \Delta\eta and \Delta\M files.
## Please note the current version of the software expects these files to NOT include the reference codon
## If you files contain the reference codon, please remove them
## Simulations can take a bit of time and memory...

library(AnaCoDa)


sel  <- "../Results/mp_sp_pseudo_sanity_check/chain/run_15/Parameter_est/sp_main_Selection"
mut <-  "../Results/mp_sp_pseudo_sanity_check/chain/run_15/Parameter_est/sp_main_Mutation"

genome.sp <- initializeGenomeObject("../Data/Genomes/Signal_peptides/sp_main.fasta")
phi <- read.table("../Data/Empirical/ROC_Phi/Signal_peptides/sp_main_phi.csv",sep=",",header=T)
parameter.sp <- initializeParameterObject(genome.sp,sphi=c(0.01),num.mixtures = 1,gene.assignment = rep(1,length(genome.sp)),mixture.definition = "allUnique",initial.expression.values =phi[,2],split.serine = TRUE)
parameter.sp$initSelectionCategories(sel,1)
parameter.sp$initMutationCategories(mut,1)
model.sp <- initializeModelObject(parameter.sp,model="ROC",with.phi=FALSE)

genome.nosp <- initializeGenomeObject("../Data/Genomes/nosp_main_first_23.fasta")
phi <- read.table("../Data/Empirical/ROC_Phi/nosp_main_first_23_phi.csv",sep=",",header=T)
parameter.nosp <- initializeParameterObject(genome.nosp,sphi=c(0.01),num.mixtures = 1,gene.assignment = rep(1,length(genome.nosp)),mixture.definition = "allUnique",initial.expression.values =10^phi[,2],split.serine = TRUE)
parameter.nosp$initSelectionCategories(sel,1)
parameter.nosp$initMutationCategories(mut,1)
model.nosp <- initializeModelObject(parameter.nosp,model="ROC",with.phi=FALSE)
# 
genome.psp <- initializeGenomeObject("../Data/Genomes/Pseudo_signal/pseudo_sp_main_norm.fasta")
phi <- read.table("../Data/Empirical/ROC_Phi/Pseudo_signal/pseudo_sp_main_norm_phi.csv",sep=",",header=T)
parameter.psp <- initializeParameterObject(genome.psp,sphi=c(0.01),num.mixtures = 1,gene.assignment = rep(1,length(genome.psp)),mixture.definition = "allUnique",initial.expression.values =phi[,2],split.serine = TRUE)
parameter.psp$initSelectionCategories(sel,1)
parameter.psp$initMutationCategories(mut,1)
model.psp <- initializeModelObject(parameter.psp,model="ROC",with.phi=FALSE)

for (i in 1:500)
{
  genome.nosp <- initializeGenomeObject("../Data/Genomes/nosp_main_first_23.fasta")
  model.nosp$simulateGenome(genome.nosp)
  genome.nosp$writeFasta(paste0("../Data/Genomes/Simulated/Simulated_nosp/nosp_sim_sp_dEta",i,".fasta"),simulated=T)
  genome.sp <- initializeGenomeObject("../Data/Genomes/Signal_peptides/sp_main.fasta")
  model.sp$simulateGenome(genome.sp)
  genome.sp$writeFasta(paste0("../Data/Genomes/Simulated/Simulated_sp/sp_sim_sp_dEta_",i,".fasta"),simulated=T)
  genome.psp <- initializeGenomeObject("../Data/Genomes/Pseudo_signal/pseudo_sp_main_norm.fasta")
  model.psp$simulateGenome(genome.psp)
  genome.psp$writeFasta(paste0("../Data/Genomes/Simulated/Simulated_psp/psp_sim_sp_dEta",i,".fasta"),simulated=T)
   
}
