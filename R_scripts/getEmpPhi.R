library(plyr)
library(ggplot2)

normalize<-function(phi)
{
  mean.phi <- mean(phi,na.rm = T)
  norm.phi <- phi/mean.phi
  return(norm.phi)
}

geom_mean<-function(x)
{
  len <- length(which(x != 0 & is.na(x)==F))
  if (len != 0){
    total = 1.0
    for(i in x[which(x != 0 & is.na(x)==F)])
    {
      total <- total * i
    }
    total <- total ^ (1/len)
  } else{
    total <- 0.0
  }
  return(total)
}

li_mrna <- read.table("../Data/Empirical/Li/li_mrna_level_main_ht_liberal.csv",sep=",",header=T,stringsAsFactors = F)
li_au <- read.table("../Data/Empirical/Li/trans_eff_main_ht_liberal.csv",sep=",",header=T,stringsAsFactors = F)
li_mops_complete <- read.table("../Data/Empirical/Li/li_prot_synth_main_ht_liberal_mops_complete.csv",sep=",",header=T,stringsAsFactors = F)
li_mops_min <- read.table("../Data/Empirical/Li/li_prot_synth_main_ht_liberal_mops_minimal.csv",sep=",",header=T,stringsAsFactors = F)
li_mops_no_meth <- read.table("../Data/Empirical/Li/li_prot_synth_main_ht_liberal_mops_complete_no_meth.csv",sep=",",header=T,stringsAsFactors = F)
cho <- read.table("../Data/Empirical/Cho/cho_mrna_level_main_ht_liberal.csv",sep=",",header=T,stringsAsFactors = F)
li <- join(li_mops_complete,li_mops_min,by="Gene")
li <- join(li,li_mops_no_meth,by="Gene")
df <- data.frame(li_mrna[,2]*li_au[,2],li_mops_complete[,2],li_mops_min[,2],li_mops_no_meth[,2],row.names = li[,1],stringsAsFactors = F)

for (i in colnames(df))
{

  df[,i] <- normalize(df[,i])
}

x<-apply(df,geom_mean,MARGIN = 1)
x[x==0]<-0.0001
x <- normalize(x)
df<-cbind(df,x)
write.table(x=df[,c(5)],file="../Data/Empirical/Emp_Phi/Ecoli_emp_phi.csv",quote=F,row.names = rownames(df), col.names = c("Phi"),sep=",")
write.table(df[1:3358,c(5)],file="../Data/Empirical/Emp_Phi/Ecoli_main_emp_phi.csv",quote = F,row.names = rownames(df)[1:3358] ,col.names = c("Phi"),sep=",")
write.table(df[3359:4140,c(5)],file="../Data/Empirical/Emp_Phi/Ecoli_pht_emp_phi.csv",quote = F,row.names =rownames(df)[3359:4140] ,col.names = c("Phi"),sep=",")
