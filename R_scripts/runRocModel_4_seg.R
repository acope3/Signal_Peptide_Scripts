library(AnaCoDa)
rm(list=ls())




#################################### Calculating Effective Sample Sizes ##############################################################
calcEffectiveSampleSize <- function(parameter,paramType,burn.in=1000,thin=100,mixture=1)
{
  trace <- parameter$getTraceObject()
  names.aa<- aminoAcids()
  for (aa in names.aa)
  {
    if (aa == "M" || aa == "W" || aa == "X")
      next
    codons <- AAToCodon(aa,TRUE)
    for (i in 1:length(codons))
    {
      for (j in 1:mixture)
      {
        csp.trace <-trace$getCodonSpecificParameterTraceByMixtureElementForCodon(j,codons[i],paramType,TRUE)
        csp.trace <- csp.trace[burn.in:length(csp.trace)]
        y <- mcmc(csp.trace,thin)
        cat(aa,codons[i],effectiveSize(y),"Mixture",j,"\n",sep = "\t")
      }
    }
  }
}
######################################################################################################################################






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
  xseg <- "mp_main"
  yseg <- "sp_main"
  zseg <- "pseudo_mp_main_test"
  wseg <- "pseudo_sp_main_test"
  mix.name.1 <- "Mature Peptide" 
  mix.name.2 <- "Signal Peptide"
  mix.name.3 <- "Pseudo Mp"
  mix.name.4 <- "Pseudo Sp"
}else{
  for(i in 1:length(args))
  {
    eval(parse(text=args[[i]]))
  }
}

genome.file.1 <- paste0("../../Data/Genomes/Ecoli_K12_ncbi/main_ht/",xseg,".fasta")
genome.file.2 <- paste0("../../Data/Genomes/Ecoli_K12_ncbi/main_ht/",yseg,".fasta")
genome.file.3 <- paste0("../../Data/Genomes/Ecoli_K12_ncbi/main_ht/skewnorm/",zseg,".fasta")
genome.file.4 <- paste0("../../Data/Genomes/Ecoli_K12_ncbi/main_ht/skewnorm/",wseg,".fasta")

expression.file.1 <- paste0("../../Data/Empirical/Ecoli_K12_ncbi/main_ht/",xseg,"_expression_improved.csv")
expression.file.2 <- paste0("../../Data/Empirical/Ecoli_K12_ncbi/main_ht/",yseg,"_expression_improved.csv")
expression.file.3 <- paste0("../../Data/Empirical/Ecoli_K12_ncbi/main_ht/skewnorm/",zseg,"_expression.csv")
expression.file.4 <- paste0("../../Data/Empirical/Ecoli_K12_ncbi/main_ht/skewnorm/",wseg,"_expression.csv")

with.phi <- FALSE

if (with.phi) {
  genome.tmp <- initializeGenomeObject(file=genome.file.1,observed.expression.file=expression.file.1,match.expression.by.id = FALSE)
  size.tmp <- length(genome.tmp)
  genome.tmp <- initializeGenomeObject(file=genome.file.2,genome = genome.tmp,observed.expression.file=expression.file.2,match.expression.by.id = FALSE,append=TRUE)
  size.tmp.2 <- length(genome.tmp)
  genome.tmp <- initializeGenomeObject(file=genome.file.3,genome = genome.tmp,observed.expression.file=expression.file.3,match.expression.by.id = FALSE,append=TRUE)
  size.tmp.3 <- length(genome.tmp)
  genome<- initializeGenomeObject(file=genome.file.4,genome = genome.tmp,observed.expression.file=expression.file.4,match.expression.by.id = FALSE,append=TRUE)
  
  } else {
   genome.tmp <- initializeGenomeObject(file=genome.file.1)
  size.tmp <- length(genome.tmp)
  genome.tmp <- initializeGenomeObject(file=genome.file.2,genome = genome.tmp,append=TRUE)
  size.tmp.2 <- length(genome.tmp)
  genome.tmp <- initializeGenomeObject(file=genome.file.3,genome = genome.tmp,append=TRUE)
  size.tmp.3 <- length(genome.tmp)
  genome<- initializeGenomeObject(file=genome.file.4,genome = genome.tmp,append=TRUE)
  
}

cat("Genome loaded\n")
#initialize parameter object
numMixtures <- 4
sphi_init <- rep(0.01,numMixtures)


mixDef <- "mutationShared"
size <- length(genome)
cat(size,"\n")
index <- c(1:size)

geneAssignment <- c(rep(1,size.tmp),rep(2,size.tmp.2-size.tmp),rep(3,size.tmp.3-size.tmp.2),rep(4,size-size.tmp.3))

segment_1_exp <- read.table(file=expression.file.1,sep=",",header=TRUE)
segment_2_exp <- read.table(file=expression.file.2,sep=",",header=TRUE)
segment_3_exp <- read.table(file=expression.file.3,sep=",",header=TRUE)
segment_4_exp <- read.table(file=expression.file.4,sep=",",header=TRUE)

init_phi <- c(segment_1_exp[,2],segment_2_exp[,2],segment_3_exp[,2],segment_4_exp[,2])
#parameter <- initializeParameterObject(genome,sphi_init,numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef,initial.expression.values = init_phi)
parameter <- initializeParameterObject(genome,sphi_init,init.sepsilon = c(0.01),numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef, initial.expression.values = init_phi)
#parameter <- initializeParameterObject(init.with.restart.file = "../../Results/Ecoli_K12_ncbi/mp_sp_pseudo_test/chain/run_3/Restart_files/rstartFile.rst_final")
# parameter$initMutationCategories(c("../../Results/Ecoli_K12_ncbi/Segmented_main_ht_fixed/chainundefined/run_2/Parameter_est/mutation_mod_Ecoli_K12_MG1655_ncbi_main_liberal.csv"),1)
# parameter$fixDM()

# initialize MCMC object
samples <-samp
thinning <- thin
adaptiveWidth <-adapt
mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth, 
                             est.expression=FALSE, est.csp=TRUE, est.hyper=FALSE,est.mix = FALSE)
mcmc$setStepsToAdapt((samples*thinning)/2)
# get model object
model <- initializeModelObject(parameter, "ROC", with.phi,fix.observation.noise = TRUE)

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
  runMCMC(mcmc, genome, model, 2,divergence.iteration = div)
)

writeParameterObject(parameter,paste(dir_name,"R_objects/parameter.Rda",sep="/"))
writeMCMCObject(mcmc,file=paste(dir_name,"R_objects/mcmc.Rda",sep="/"))


sel.1 <-paste0("Parameter_est/selection_",xseg,".csv")
sel.2 <-paste0("Parameter_est/selection_",yseg,".csv")
sel.3 <-paste0("Parameter_est/selection_",zseg,".csv")
sel.4 <-paste0("Parameter_est/selection_",wseg,".csv")
mut.1 <-paste0("Parameter_est/mutation_",xseg,".csv")
mut.2 <-paste0("Parameter_est/mutation_",yseg,".csv")
mut.3 <-paste0("Parameter_est/mutation_",zseg,".csv")
mut.4 <-paste0("Parameter_est/mutation_",wseg,".csv")

getCSPEstimates(parameter,paste(dir_name,sel.1,sep="/"),"Selection",1,samples*0.5)
getCSPEstimates(parameter,paste(dir_name,mut.1,sep="/"),"Mutation",1,samples*0.5)
getCSPEstimates(parameter,paste(dir_name,sel.2,sep="/"),"Selection",2,samples*0.5)
getCSPEstimates(parameter,paste(dir_name,mut.2,sep="/"),"Mutation",2,samples*0.5)
getCSPEstimates(parameter,paste(dir_name,sel.3,sep="/"),"Selection",3,samples*0.5)
getCSPEstimates(parameter,paste(dir_name,mut.3,sep="/"),"Mutation",3,samples*0.5)
getCSPEstimates(parameter,paste(dir_name,sel.4,sep="/"),"Selection",4,samples*0.5)
getCSPEstimates(parameter,paste(dir_name,mut.4,sep="/"),"Mutation",4,samples*0.5)

# mixtureAssignment <- getMixtureAssignmentEstimate(parameter,c(1:size),samples*0.25)
# expressionValues <- getExpressionEstimatesForMixture(parameter,c(1:size),mixtureAssignment,samples*0.25)
# expressionValues <- log10(expressionValues)
# write(expressionValues,file=paste(dir_name,"Parameter_est/gene_expression.txt",sep="/"),ncolumns = 1)


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
plot(parameter,what="Mutation",samples=samples*0.5,mixture.name=c(mix.name.1,mix.name.2,mix.name.3,mix.name.4))
plot(parameter,what="Selection",samples=samples*0.5,mixture.name=c(mix.name.1,mix.name.2,mix.name.3,mix.name.4))
plot(trace, what = "Mutation", mixture = 1)
plot(trace, what = "Selection", mixture = 1)
plot(trace, what = "Mutation", mixture = 2)
plot(trace, what = "Selection", mixture = 2)
plot(trace, what = "Mutation", mixture = 3)
plot(trace, what = "Selection", mixture = 3)
plot(trace, what = "Mutation", mixture = 4)
plot(trace, what = "Selection", mixture = 4)
plot(model, genome, samples = samples*0.5, mixture = 1,main = mix.name.1)
plot(model, genome, samples = samples*0.5, mixture = 2,main = mix.name.2)
plot(model, genome, samples = samples*0.5, mixture = 3,main = mix.name.3)
plot(model, genome, samples = samples*0.5, mixture = 4,main = mix.name.4)
dev.off()


diag <- convergence.test(mcmc,samples = samples*0.5,thin=thinning)
z<-abs(diag$z)
while(z>1.96 || done.adapt==FALSE)
{
  samples <- 5000
  parameter<-initializeParameterObject(init.with.restart.file = paste(dir_name,"Restart_files/rstartFile.rst_final",sep="/"))
  #parameter$fixDM()
  run_number <- run_number + 1
  dir_name <- paste0(directory,"/run_",run_number)
  dir.create(dir_name)
  dir.create(paste(dir_name,"Graphs",sep="/"))
  dir.create(paste(dir_name,"Restart_files",sep="/"))
  dir.create(paste(dir_name,"Parameter_est",sep="/"))
  dir.create(paste(dir_name,"R_objects",sep="/"))
  
  mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth, 
                               est.expression=FALSE, est.csp=TRUE, est.hyper=FALSE,est.mix=FALSE)
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
    runMCMC(mcmc, genome, model, 2)
  )
  
  writeParameterObject(parameter,paste(dir_name,"R_objects/parameter.Rda",sep="/"))
  writeMCMCObject(mcmc,file=paste(dir_name,"R_objects/mcmc.Rda",sep="/"))
  
  getCSPEstimates(parameter,paste(dir_name,sel.1,sep="/"),"Selection",1,samples*0.5)
  getCSPEstimates(parameter,paste(dir_name,mut.1,sep="/"),"Mutation",1,samples*0.5)
  getCSPEstimates(parameter,paste(dir_name,sel.2,sep="/"),"Selection",2,samples*0.5)
  getCSPEstimates(parameter,paste(dir_name,mut.2,sep="/"),"Mutation",2,samples*0.5)
  getCSPEstimates(parameter,paste(dir_name,sel.3,sep="/"),"Selection",3,samples*0.5)
  getCSPEstimates(parameter,paste(dir_name,mut.3,sep="/"),"Mutation",3,samples*0.5)
  getCSPEstimates(parameter,paste(dir_name,sel.4,sep="/"),"Selection",4,samples*0.5)
  getCSPEstimates(parameter,paste(dir_name,mut.4,sep="/"),"Mutation",4,samples*0.5)
  
  # mixtureAssignment <- getMixtureAssignmentEstimate(parameter,c(1:size),samples*0.25)
  # expressionValues <- getExpressionEstimatesForMixture(parameter,c(1:size),mixtureAssignment,samples*0.25)
  # expressionValues <- log10(expressionValues)
  # write(expressionValues,file=paste(dir_name,"Parameter_est/gene_expression.txt",sep="/"),ncolumns = 1)
  
  
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
  plot(parameter,what="Mutation",samples=samples*0.5,mixture.name=c(mix.name.1,mix.name.2,mix.name.3,mix.name.4))
  plot(parameter,what="Selection",samples=samples*0.5,mixture.name=c(mix.name.1,mix.name.2,mix.name.3,mix.name.4))
  plot(trace, what = "Mutation", mixture = 1)
  plot(trace, what = "Selection", mixture = 1)
  plot(trace, what = "Mutation", mixture = 2)
  plot(trace, what = "Selection", mixture = 2)
  plot(trace, what = "Mutation", mixture = 3)
  plot(trace, what = "Selection", mixture = 3)
  plot(trace, what = "Mutation", mixture = 4)
  plot(trace, what = "Selection", mixture = 4)
  plot(model, genome, samples = samples*0.5, mixture = 1,main = mix.name.1)
  plot(model, genome, samples = samples*0.5, mixture = 2,main = mix.name.2)
  plot(model, genome, samples = samples*0.5, mixture = 3,main = mix.name.3)
  plot(model, genome, samples = samples*0.5, mixture = 4,main = mix.name.4)
  dev.off()
  
  
  diag <- convergence.test(mcmc,samples = samples*0.5,thin=thinning)
  z<-abs(diag$z)
}

samples <- samp
thinning <- thin
parameter<-initializeParameterObject(init.with.restart.file = paste(dir_name,"Restart_files/rstartFile.rst_final",sep="/"))
#parameter$fixDM()
run_number <- run_number + 1
dir_name <- paste0(directory,"/run_",run_number)
dir.create(dir_name)
dir.create(paste(dir_name,"Graphs",sep="/"))
dir.create(paste(dir_name,"Restart_files",sep="/"))
dir.create(paste(dir_name,"Parameter_est",sep="/"))
dir.create(paste(dir_name,"R_objects",sep="/"))

mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth, 
                             est.expression=FALSE, est.csp=TRUE, est.hyper=FALSE,est.mix=FALSE)

mcmc$setStepsToAdapt(0)

model <- initializeModelObject(parameter, "ROC", with.phi)
setRestartSettings(mcmc, paste(dir_name,"Restart_files/rstartFile.rst",sep="/"), adaptiveWidth, TRUE)
#run mcmc on genome with parameter using model
system.time(
  runMCMC(mcmc, genome, model, 2)
)

writeParameterObject(parameter,paste(dir_name,"R_objects/parameter.Rda",sep="/"))
writeMCMCObject(mcmc,file=paste(dir_name,"R_objects/mcmc.Rda",sep="/"))


getCSPEstimates(parameter,paste(dir_name,sel.1,sep="/"),"Selection",1,samples)
getCSPEstimates(parameter,paste(dir_name,mut.1,sep="/"),"Mutation",1,samples)
getCSPEstimates(parameter,paste(dir_name,sel.2,sep="/"),"Selection",2,samples)
getCSPEstimates(parameter,paste(dir_name,mut.2,sep="/"),"Mutation",2,samples)
getCSPEstimates(parameter,paste(dir_name,sel.3,sep="/"),"Selection",3,samples)
getCSPEstimates(parameter,paste(dir_name,mut.3,sep="/"),"Mutation",3,samples)
getCSPEstimates(parameter,paste(dir_name,sel.4,sep="/"),"Selection",4,samples)
getCSPEstimates(parameter,paste(dir_name,mut.4,sep="/"),"Mutation",4,samples)

# mixtureAssignment <- getMixtureAssignmentEstimate(parameter,c(1:size),samples*0.25)
# expressionValues <- getExpressionEstimatesForMixture(parameter,c(1:size),mixtureAssignment,samples*0.25)
# expressionValues <- log10(expressionValues)
# write(expressionValues,file=paste(dir_name,"Parameter_est/gene_expression.txt",sep="/"),ncolumns = 1)


#plots different aspects of trace
trace <- parameter$getTraceObject()
pdf(paste(dir_name,"Graphs/mcmc_traces.pdf",sep="/"))
plot(mcmc,what = "LogPosterior")
plot(trace, what = "ExpectedPhi")
dev.off()

pdf(paste(dir_name,"Graphs/CSP_traces_CUB_plot.pdf",sep="/"), width = 11, height = 12)
plot(parameter,what="Mutation",samples=samples,mixture.name=c(mix.name.1,mix.name.2,mix.name.3,mix.name.4))
plot(parameter,what="Selection",samples=samples,mixture.name=c(mix.name.1,mix.name.2,mix.name.3,mix.name.4))
plot(trace, what = "Mutation", mixture = 1)
plot(trace, what = "Selection", mixture = 1)
plot(trace, what = "Mutation", mixture = 2)
plot(trace, what = "Selection", mixture = 2)
plot(trace, what = "Mutation", mixture = 3)
plot(trace, what = "Selection", mixture = 3)
plot(trace, what = "Mutation", mixture = 4)
plot(trace, what = "Selection", mixture = 4)
plot(model, genome, samples = samples, mixture = 1,main = mix.name.1)
plot(model, genome, samples = samples, mixture = 2,main = mix.name.2)
plot(model, genome, samples = samples, mixture = 3,main = mix.name.3)
plot(model, genome, samples = samples, mixture = 4,main = mix.name.4)
dev.off()
