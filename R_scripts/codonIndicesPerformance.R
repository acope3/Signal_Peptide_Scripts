## Author: Alex Cope
## Script for analyzing CUB
## Contains a modified form of the CAI function found in AnaCoDa that should reduce potential numerical issues.
## Some filepaths are hardcoded and will need to be changed.
## The file ecolik12.trna can be found on the github page for the tAI package

library(AnaCoDa)
library(tAI)

getCAIweights <- function(referenceGenome)
{
  aa.vec <- aminoAcids()
  aa.vec <- aa.vec[-length(aa.vec)]
  
  wi.list <- vector(mode = "list", length = length(aa.vec))
  names(wi.list) <- aa.vec
  
  codon.names <- NULL
  for(aa in aa.vec)
  {
    codon.names <- c(codon.names, AAToCodon(aa))
    ## create reference table for each codon and gene
    codonCountForAA.ref <- AnaCoDa:::getCodonCountsForAA(aa, genome = referenceGenome)
    fi <- colSums(codonCountForAA.ref)
    fi.max <- max(fi)
    wi.list[[aa]] <- fi / fi.max
  }
  
  wi.vec <- unlist(wi.list)
  wi.vec[wi.vec == 0.0] <- 0.0001
  names(wi.vec) <- codon.names
  return(wi.vec)
}


### NOT EXPOSED
calcCAI <- function(gene, wi)
{
  # sequence string to triplets
  seq <- gene$seq
  seq <- unlist(strsplit(seq, ""))
  seq <- paste(seq[c(T,F,F)], seq[c(F,T,F)], seq[c(F,F,T)], sep="")
  codon.length <- length(seq)
  
  CAI <- 0
  for(s in seq)
  {
    if(is.na(wi[s]) || s == "ATG" || s == "TGG" || s == "TAG" || s=="TAA" || s=="TGA") 
      {
        codon.length <- codon.length - 1
        next
    }    
      CAI <- CAI + log(wi[s])
  }
  CAI <- exp((1/codon.length)*CAI)
  return(CAI)
}


getCAI <- function(referenceGenome, testGenome)
{
  genes <- testGenome$getGenes(FALSE)
  wi <- getCAIweights(referenceGenome)
  CAI <- unlist(lapply(genes, calcCAI, wi))
  names(CAI) <- getNames(testGenome, FALSE)
  return(CAI)  
}



ribo<- initializeGenomeObject("../Data/Genomes/ribo_prot.fasta")
pval<-c()
sp.total <- c()
nosp.total <- c()
for (i in 1:500)
{
  cat(i,"\n")
  sp <- initializeGenomeObject(paste0("../Data/Genomes/Simulated/Simulated_sp/sp_sim_sp_dEta_",i,".fasta"))
  nosp <- initializeGenomeObject(paste0("../Data/Genomes/Simulated/Simulated_psp/psp_sim_sp_dEta",i,".fasta"))
  sp.cai <- getCAI(ribo,sp)
  nosp.cai <- getCAI(ribo,nosp)
  sp.total <- c(sp.total,sp.cai)
  nosp.total <- c(nosp.total,nosp.cai)
  x<-t.test(sp.cai,nosp.cai,alternative="less")
  pval<-c(pval,x$p.value)
}
hist(pval,breaks=seq(0,1,0.05),ylim = c(0,500))
write(pval,file = "../pval.txt",ncolumns = 1)
pval<-c()
sp.total <- c()
nosp.total <- c()
eco.trna <- scan("~/tai/misc/ecolik12.trna")
eco.ws <- get.ws(tRNA = eco.trna,sking=1)
for (i in 1:500)
{
  cat(i,"\n")
  sp.m <- matrix(scan(paste0("../Data/Genomes/Simulated/Simulated_sp/sp_sim_",i,".m")),ncol=61,byrow=TRUE)
  nosp.m <- matrix(scan(paste0("../Data/Genomes/Simulated/Simulated_nosp_aa_norm/nosp_sim_",i,".m")),ncol=61,byrow=TRUE)
  sp.m <- sp.m[,-33]
  nosp.m <- nosp.m[,-33]
  sp.tai <- get.tai(sp.m,eco.ws)
  nosp.tai <- get.tai(nosp.m,eco.ws)
  sp.total <- c(sp.total,sp.tai)
  nosp.total <- c(nosp.total,nosp.tai)
  x<-t.test(sp.tai,nosp.tai,alternative="greater")
  pval<-c(pval,x$p.value)
}
hist(pval,breaks=seq(0,1,0.05))
write(pval,file = "../pval_tai.txt",ncolumns = 1)


