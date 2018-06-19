library(plyr)

df <- read.table("samples.csv",sep=",",header=T)
phi <- read.table("../main_ht/gene_expression_liberal.txt",sep="",header=F)

pp <- c()
ps <- c()
pdf("all_possible_comp_mrna_abundance.pdf")
for (i in 2:(ncol(df)-1))
{
  for (j in (i+1):ncol(df))
  {
    x <- which(df[,i] != 0)
    y <- which(df[,j] !=0 )
    z <- intersect(x,y)
    plot(log10(df[z,i]),log10(df[z,j]),main=bquote("Sample " ~ .(i-1) ~" vs Sample " ~ .(j-1)),xlab=bquote("Sample " ~ .(i-1) ~ ", log10(mRNA abundance)"),ylab=bquote("Sample " ~ .(j-1) ~ ", log10(mRNA abundance)"))
    rho.p <- cor(log10(df[z,i]),log10(df[z,j]))
    rho.s <- cor(log10(df[z,i]),log10(df[z,j]),method="spearman")
    xlim <- c(min(log10(df[z,i])),max(log10(df[z,i])))
    ylim <- c(min(log10(df[z,j])),max(log10(df[z,j])))
    text(xlim[1]+0.5,ylim[2]-0.25,labels=bquote(rho["p"] ~ "=" ~ .(rho.p)))
  #  text(xlim[1]+0.5,ylim[2]-0.5,labels=bquote(rho["s"] ~ "=" ~ .(rho.s)))
    pp <- c(pp,rho.p)
    ps <- c(ps,rho.s)
  }
}
hist(x = pp,nclass = 15,main = "Pearson Correlation Coefficients\nAccession GSE48829",xlab="Pearson correlation coefficient")
#hist(x = ps,nclass = 15,main = "Spearman Correlation Coefficients\nAccession GSE48829",xlab="Spearman correlation coefficient")
dev.off()