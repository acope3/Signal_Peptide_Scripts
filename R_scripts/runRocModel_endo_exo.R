## Author: Alex Cope
## Script for analyzing CUB with endogenous and exogenous genes treated as separate categories
## Produces \phi values
## Note some function calls in this script may be out of date with the current version of AnaCoDa

library(AnaCoDa)
rm(list=ls())


CHAIN_ID <- Sys.getenv("SGE_TASK_ID")
cat("Chain ID is ",CHAIN_ID)
args<-(commandArgs(TRUE))

if(length(args)==0)
{
  div <- 0 
  directory <- getwd()
  thin <- 100
  adapt <- 100
  samp <- 10000
  xseg <- "signal_peptide_nt"
  yseg <- "mature_peptide_nt"
  mix.name.1 <- "Signal Peptide" 
  mix.name.2 <- "Mature Peptide"
}else{
  for(i in 1:length(args))
  {
    eval(parse(text=args[[i]]))
  }
}

genome.file.1 <- paste0("../../Data/Genomes/Ecoli_K12_ncbi/main_ht/",xseg,".fasta")
genome.file.2 <- paste0("../../Data/Genomes/Ecoli_K12_ncbi/main_ht/",yseg,".fasta")
expression.file.1 <- paste0("../../Data/Empirical/Ecoli_K12_ncbi/main_ht/",xseg,"_emp_001.csv")
expression.file.2 <- paste0("../../Data/Empirical/Ecoli_K12_ncbi/main_ht/",yseg,"_emp_001.csv")
with.phi <- TRUE

if (with.phi) {
  genome.tmp <- initializeGenomeObject(file=genome.file.1,observed.expression.file=expression.file.1,match.expression.by.id = FALSE)
  size.tmp <- length(genome.tmp)
  genome <- initializeGenomeObject(file=genome.file.2,genome=genome.tmp,observed.expression.file=expression.file.2,append=T,match.expression.by.id = FALSE)
} else {
  genome.tmp <- initializeGenomeObject(file=genome.file.1)
  size.tmp <- length(genome.tmp)
  genome <- initializeGenomeObject(file=genome.file.2,genome = genome.tmp)
}

cat("Genome loaded\n")
#initialize parameter object
numMixtures <- 2
sphi_init <- rep(1,numMixtures)


mixDef <- "mutationShared"
size <- length(genome)
cat(size,"\n")
index <- c(1:size)

geneAssignment <- c(rep(1,size.tmp),rep(2,size-size.tmp))

s.epsilon <- rep(7,1)

parameter <- initializeParameterObject(genome,sphi_init,init.sepsilon = s.epsilon, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)
#parameter <- initializeParameterObject(genome,sphi_init,init.sepsilon = c(0.01),numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef, initial.expression.values = init_phi)
#parameter <- initializeParameterObject(init.with.restart.file = "../../Results/Ecoli_K12_ncbi/Full_genome_emp_allUnique/chain/run_2/Restart_files/rstartFile.rst_final")


# initialize MCMC object
samples <-samp
thinning <- thin
adaptiveWidth <-adapt
mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE,est.mix = FALSE)
mcmc$setStepsToAdapt((samples*thinning)/2)
# get model object
model <- initializeModelObject(parameter, "ROC", with.phi)

run_number <- 1
dir.create(directory)
if (CHAIN_ID == 1)
{
  directory <- paste0(directory,"/nodiv")
  dir.create(directory)
}else{
  directory <- paste0(directory,"/chain",CHAIN_ID)
  dir.create(directory)
}
dir_name <- paste0(directory,"/run_",run_number)
dir.create(dir_name)
dir.create(paste(dir_name,"Graphs",sep="/"))
dir.create(paste(dir_name,"Restart_files",sep="/"))
dir.create(paste(dir_name,"Parameter_est",sep="/"))
dir.create(paste(dir_name,"R_objects",sep="/"))
setRestartSettings(mcmc, paste(dir_name,"Restart_files/rstartFile.rst",sep="/"), adaptiveWidth, TRUE)
#run mcmc on genome with parameter using model
system.time(
  runMCMC(mcmc, genome, model, 8,divergence.iteration = div)
)

writeParameterObject(parameter,paste(dir_name,"R_objects/parameter.Rda",sep="/"))
writeMCMCObject(mcmc,file=paste(dir_name,"R_objects/mcmc.Rda",sep="/"))


sel.1 <-paste0("Parameter_est/selection_",xseg,".csv")
sel.2 <-paste0("Parameter_est/selection_",yseg,".csv")
mut.1 <-paste0("Parameter_est/mutation_",xseg,".csv")
mut.2 <-paste0("Parameter_est/mutation_",yseg,".csv")

getCSPEstimates(parameter,paste(dir_name,sel.1,sep="/"),"Selection",1,samples*0.5)
getCSPEstimates(parameter,paste(dir_name,mut.1,sep="/"),"Mutation",1,samples*0.5)
getCSPEstimates(parameter,paste(dir_name,sel.2,sep="/"),"Selection",2,samples*0.5)
getCSPEstimates(parameter,paste(dir_name,mut.2,sep="/"),"Mutation",2,samples*0.5)


# mixtureAssignment <- getMixtureAssignmentEstimate(parameter,c(1:size),samples*0.25)
# expressionValues <- getExpressionEstimates(parameter,c(1:size),samples*0.25)
# write.table(expressionValues,file=paste(dir_name,"Parameter_est/gene_expression.txt",sep="/"),quote = F,sep=",")
mixtureAssignment <- getMixtureAssignmentEstimate(parameter,c(1:size),samples*0.5)
expressionValues <- getExpressionEstimates(parameter,c(1:size),samples*0.5)
write.table(expressionValues,file=paste(dir_name,"Parameter_est/gene_expression.txt",sep="/"),sep=",",col.names = T,quote = F,row.names = F)


#plots different aspects of trace
trace <- parameter$getTraceObject()
pdf(paste(dir_name,"Graphs/mcmc_traces.pdf",sep="/"))
plot(mcmc,what = "LogPosterior")
plot(trace, what = "ExpectedPhi")
aa <- aminoAcids()
done.adapt <- TRUE
for(a in aa)
{
  if (a=="M"||a=="X"||a=="W") next
  accept.trace <- trace$getCodonSpecificAcceptanceRateTraceForAA(a)
  len <- length(accept.trace)
  mean.acceptance <- mean(accept.trace[(len-len*0.1):len])
  if (mean.acceptance < 0.1 || mean.acceptance > 0.44) done.adapt <- FALSE
  plot(accept.trace,main=paste0("Acceptace Rate for ",a),xlab="Samples",ylab="Acceptance Rate",type="l")
}
dev.off()

pdf(paste(dir_name,"Graphs/CSP_traces_CUB_plot.pdf",sep="/"), width = 11, height = 12)
plot(parameter,what="Mutation",samples=samples*0.5,mixture.name=c(mix.name.1,mix.name.2))
plot(parameter,what="Selection",samples=samples*0.5,mixture.name=c(mix.name.1,mix.name.2))
plot(trace, what = "Mutation", mixture = 1)
plot(trace, what = "Selection", mixture = 1)
plot(trace, what = "Mutation", mixture = 2)
plot(trace, what = "Selection", mixture = 2)
plot(model, genome, samples = samples*0.5, mixture = 1,main = mix.name.1)
plot(model, genome, samples = samples*0.5, mixture = 2,main = mix.name.2)
dev.off()


diag <- convergence.test(mcmc,samples = samples*0.5,thin=thinning)
z<-abs(diag$z)
while(z>1.96 || done.adapt==FALSE)
{
  samples <- 5000
  parameter<-initializeParameterObject(init.with.restart.file = paste(dir_name,"Restart_files/rstartFile.rst_final",sep="/"))
  run_number <- run_number + 1
  dir_name <- paste0(directory,"/run_",run_number)
  dir.create(dir_name)
  dir.create(paste(dir_name,"Graphs",sep="/"))
  dir.create(paste(dir_name,"Restart_files",sep="/"))
  dir.create(paste(dir_name,"Parameter_est",sep="/"))
  dir.create(paste(dir_name,"R_objects",sep="/"))
  
  mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth, 
                               est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE,est.mix=FALSE)
  if(!done.adapt)
  {
    mcmc$setStepsToAdapt((samples*thinning)/2)
  } else{
    mcmc$setStepsToAdapt(0)
  }
  model <- initializeModelObject(parameter, "ROC", with.phi)
  setRestartSettings(mcmc, paste(dir_name,"Restart_files/rstartFile.rst",sep="/"), adaptiveWidth, TRUE)
  #run mcmc on genome with parameter using model
  system.time(
    runMCMC(mcmc, genome, model, 8)
  )
  
  writeParameterObject(parameter,paste(dir_name,"R_objects/parameter.Rda",sep="/"))
  writeMCMCObject(mcmc,file=paste(dir_name,"R_objects/mcmc.Rda",sep="/"))
  
  getCSPEstimates(parameter,paste(dir_name,sel.1,sep="/"),"Selection",1,samples*0.5)
  getCSPEstimates(parameter,paste(dir_name,mut.1,sep="/"),"Mutation",1,samples*0.5)
  getCSPEstimates(parameter,paste(dir_name,sel.2,sep="/"),"Selection",2,samples*0.5)
  getCSPEstimates(parameter,paste(dir_name,mut.2,sep="/"),"Mutation",2,samples*0.5)
  
  mixtureAssignment <- getMixtureAssignmentEstimate(parameter,c(1:size),samples*0.5)
  expressionValues <- getExpressionEstimates(parameter,c(1:size),samples*0.5)
  write.table(expressionValues,file=paste(dir_name,"Parameter_est/gene_expression.txt",sep="/"),sep=",",col.names = T,quote = F,row.names = F)
  
  trace <- parameter$getTraceObject()
  pdf(paste(dir_name,"Graphs/mcmc_traces.pdf",sep="/"))
  plot(mcmc,what = "LogPosterior")
  plot(trace, what = "ExpectedPhi")
  if(!done.adapt)
  {
    aa <- aminoAcids()
    done.adapt <- TRUE
    for(a in aa)
    {
      if (a=="M"||a=="X"||a=="W") next
      accept.trace <- trace$getCodonSpecificAcceptanceRateTraceForAA(a)
      len <- length(accept.trace)
      mean.acceptance <- mean(accept.trace[(len-len*0.1):len])
      if (mean.acceptance < 0.1 || mean.acceptance > 0.44) done.adapt <- FALSE
      plot(accept.trace,main=paste0("Acceptace Rate for ",a),xlab="Samples",ylab="Acceptance Rate",type="l")
    }
  }
  dev.off()
  
  pdf(paste(dir_name,"Graphs/CSP_CUB_plot.pdf",sep="/"), width = 11, height = 12)
  plot(parameter,what="Mutation",samples=samples*0.5,mixture.name=c(mix.name.1,mix.name.2))
  plot(parameter,what="Selection",samples=samples*0.5,mixture.name=c(mix.name.1,mix.name.2))
  plot(trace, what = "Mutation", mixture = 1)
  plot(trace, what = "Selection", mixture = 1)
  plot(trace, what = "Mutation", mixture = 2)
  plot(trace, what = "Selection", mixture = 2)
  plot(model, genome, samples = samples*0.5, mixture = 1,main = mix.name.1)
  plot(model, genome, samples = samples*0.5, mixture = 2,main = mix.name.2)
  dev.off()
  
  
  diag <- convergence.test(mcmc,samples = samples*0.5,thin=thinning)
  z<-abs(diag$z)
}

samples <- 5000
thinning <- thin
parameter<-initializeParameterObject(init.with.restart.file = paste(dir_name,"Restart_files/rstartFile.rst_final",sep="/"))
run_number <- run_number + 1
dir_name <- paste0(directory,"/run_",run_number)
dir.create(dir_name)
dir.create(paste(dir_name,"Graphs",sep="/"))
dir.create(paste(dir_name,"Restart_files",sep="/"))
dir.create(paste(dir_name,"Parameter_est",sep="/"))
dir.create(paste(dir_name,"R_objects",sep="/"))

mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE,est.mix=FALSE)

mcmc$setStepsToAdapt(0)

model <- initializeModelObject(parameter, "ROC", with.phi)
setRestartSettings(mcmc, paste(dir_name,"Restart_files/rstartFile.rst",sep="/"), adaptiveWidth, TRUE)
#run mcmc on genome with parameter using model
system.time(
  runMCMC(mcmc, genome, model, 8)
)

writeParameterObject(parameter,paste(dir_name,"R_objects/parameter.Rda",sep="/"))
writeMCMCObject(mcmc,file=paste(dir_name,"R_objects/mcmc.Rda",sep="/"))

getCSPEstimates(parameter,paste(dir_name,sel.1,sep="/"),"Selection",1,samples)
getCSPEstimates(parameter,paste(dir_name,mut.1,sep="/"),"Mutation",1,samples)
getCSPEstimates(parameter,paste(dir_name,sel.2,sep="/"),"Selection",2,samples)
getCSPEstimates(parameter,paste(dir_name,mut.2,sep="/"),"Mutation",2,samples)

mixtureAssignment <- getMixtureAssignmentEstimate(parameter,c(1:size),samples)
expressionValues <- getExpressionEstimates(parameter,c(1:size),samples)
write.table(expressionValues,file=paste(dir_name,"Parameter_est/gene_expression.txt",sep="/"),sep=",",col.names = T,quote = F,row.names = F)


trace <- parameter$getTraceObject()
pdf(paste(dir_name,"Graphs/mcmc_traces.pdf",sep="/"))
plot(mcmc,what = "LogPosterior")
plot(trace, what = "ExpectedPhi")
dev.off()

pdf(paste(dir_name,"Graphs/CSP_CUB_plot.pdf",sep="/"), width = 11, height = 12)
plot(parameter,what="Mutation",samples=samples,mixture.name=c(mix.name.1,mix.name.2))
plot(parameter,what="Selection",samples=samples,mixture.name=c(mix.name.1,mix.name.2))
plot(trace, what = "Mutation", mixture = 1)
plot(trace, what = "Selection", mixture = 1)
plot(trace, what = "Mutation", mixture = 2)
plot(trace, what = "Selection", mixture = 2)
plot(model, genome, samples = samples, mixture = 1,main = mix.name.1)
plot(model, genome, samples = samples, mixture = 2,main = mix.name.2)
dev.off()
