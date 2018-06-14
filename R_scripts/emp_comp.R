## Author: Alex Cope
## Generate ROC \phi vs. Empirical \phi comparisons.
## Will need to update file paths!

library(plyr)
library(ggplot2)

li_mrna <- read.table("../Data/Empirical/Li/li_mrna_level_main_ht_liberal.csv",sep=",",header=T)
li_au <- read.table("../Data/Empirical/Li/trans_eff_main_ht_liberal.csv",sep=",",header=T)
li_mops_complete <- read.table("../Data/Empirical/Li/li_prot_synth_main_ht_liberal_mops_complete.csv",sep=",",header=T)
li_mops_min <- read.table("../Data/Empirical/Li/li_prot_synth_main_ht_liberal_mops_minimal.csv",sep=",",header=T)
li_mops_no_meth <- read.table("../Data/Empirical/Li/li_prot_synth_main_ht_liberal_mops_complete_no_meth.csv",sep=",",header=T)
li <- join(li_mops_complete,li_mops_min,by="Gene")
li <- join(li,li_mops_no_meth,by="Gene")
cho <- read.table("../Data/Empirical/Cho/cho_mrna_level_main_ht_liberal.csv",sep=",",header=T)
colombos <- read.table("../Data/Empirical/COLOMBOS/samples.csv",sep=",",header=T)
maksym_1<- read.table("../Data/Empirical/Maksym/samples_1_6.csv",sep=",",header=T)
maksym_2<- read.table("../Data/Empirical/Maksym/samples_7_22.csv",sep=",",header=T)
maksym <- join(maksym_1,maksym_2,by="Gene")

## Read in ROC \phi values
phi <- read.table("../Data/Empirical/ROC_Phi/Ecoli_K12_MG1655_main_pht_phi.csv",sep=",",header=T)
in.lab.cor.p <- c()
phi.cor.p <- c()


maksym.exclude <- c(1,2,3,4,5,6,7,8,17,18)
colombos.exclude <- c()
pdf("../Graphs/Maksym_in_lab_comparisons.pdf")
for (i in 2:(ncol(maksym)-1))
{
  for (j in (i+1):ncol(maksym))
  {
  
    if ((i-1) %in% maksym.exclude || j %in% maksym.exclude) next
    x <- which(maksym[1:3358,i] >= 10)
    y <- which(maksym[1:3358,j] >= 10)
    z <- intersect(x,y)
    plot(log10(maksym[z,i]),log10(maksym[z,j]),main=bquote("Sample " ~ .(i-1) ~" vs Sample " ~ .(j-1)),xlab=bquote("Sample " ~ .(i-1) ~ ", log10(mRNA abundance)"),ylab=bquote("Sample " ~ .(j-1) ~ ", log10(mRNA abundance)"))
    rho.p <- cor(log10(maksym[z,i]),log10(maksym[z,j]))
    xlim <- c(min(log10(maksym[z,i])),max(log10(maksym[z,i])))
    ylim <- c(min(log10(maksym[z,j])),max(log10(maksym[z,j])))
    text(xlim[1]+0.5,ylim[2]-0.25,labels=bquote(rho["p"] ~ "=" ~ .(rho.p)))
    in.lab.cor.p <- c(in.lab.cor.p,rho.p)
  }
}
dev.off()

pdf("../Graphs/Colombos_in_lab_comparisons.pdf")
for (i in 2:(ncol(colombos)-1))
{
  for (j in (i+1):ncol(colombos))
  {
    if ((i-1) %in% colombos.exclude || j %in% colombos.exclude) next
    x <- which(colombos[1:3358,i] >= 10)
    y <- which(colombos[1:3358,j] >= 10 )
    z <- intersect(x,y)
    plot(log10(colombos[z,i]),log10(colombos[z,j]),main=bquote("Sample " ~ .(i-1) ~" vs Sample " ~ .(j-1)),xlab=bquote("Sample " ~ .(i-1) ~ ", log10(mRNA abundance)"),ylab=bquote("Sample " ~ .(j-1) ~ ", log10(mRNA abundance)"))
    rho.p <- cor(log10(colombos[z,i]),log10(colombos[z,j]))
    xlim <- c(min(log10(colombos[z,i])),max(log10(colombos[z,i])))
    ylim <- c(min(log10(colombos[z,j])),max(log10(colombos[z,j])))
    text(xlim[1]+0.5,ylim[2]-0.25,labels=bquote(rho["p"] ~ "=" ~ .(rho.p)))
    in.lab.cor.p <- c(in.lab.cor.p,rho.p)
  }
}
dev.off()

pdf("../Graphs/Li_in_lab_comparisons.pdf")
for (i in 2:(ncol(li)-1))
{
  for (j in (i+1):ncol(li))
  {
    x <- which(li[1:3358,i] != 0)
    y <- which(li[1:3358,j] != 0)
    z <- intersect(x,y)
    if (i == 2) sample.x <- "MOPS Complete"
    else if (i == 3) sample.x <- "MOPS Minimal"
    else if (i == 4) sample.x <- "MOPS Complete, no methionine"
    if (j == 2) sample.y <- "MOPS Complete"
    else if (j == 3) sample.y <- "MOPS Minimal"
    else if (j == 4) sample.y <- "MOPS Complete, no methionine"
    plot(log10(li[z,i]),log10(li[z,j]),main=bquote(.(sample.x) ~ "vs" ~ (.sample.y)),xlab=bquote(.(sample.x) ~ ", log10(mRNA abundance)"),ylab=bquote(.(sample.y) ~ ", log10(mRNA abundance)"))
    rho.p <- cor(log10(li[z,i]),log10(li[z,j]))
    xlim <- c(min(log10(li[z,i])),max(log10(li[z,i])))
    ylim <- c(min(log10(li[z,j])),max(log10(li[z,j])))
    text(xlim[1]+0.5,ylim[2]-0.25,labels=bquote(rho["p"] ~ "=" ~ .(rho.p)))
    in.lab.cor.p <- c(in.lab.cor.p,rho.p)
  }
}


bt.lab.cor.p <- c()
pdf("../Graphs/li_cho_comp.pdf")
x <- which(li_mrna[1:3358,2]!=0)
y <- which(cho[1:3358,2] >= 10)
z <- intersect(x,y)
plot(log10(li_mrna[z,2]),log10(cho[z,2]),main="Li et al vs Cho et al:\nmRNA Abundance",xlab="Li et al, log10(mRNA abundance)",ylab="Cho et al, log10(mRNA abundance)")
rho.p <- cor(log10(li_mrna[z,2]),log10(cho[z,2]))
xlim <- c(min(log10(li_mrna[z,2])),max(log10(li_mrna[z,2])))
ylim <- c(min(log10(cho[z,2])),max(log10(cho[z,2])))
text(xlim[2]-0.5,ylim[1]+0.5,labels=bquote(rho["p"] ~ "=" ~ .(rho.p)))
bt.lab.cor.p <- c(bt.lab.cor.p,rho.p)
dev.off()

pdf("../Graphs/colombos_maksym_comp.pdf")
for (i in 2:ncol(maksym))
{
  for (j in 2:ncol(colombos))
  {
    if ((i-1) %in% maksym.exclude || (j-1) %in% colombos.exclude) next
    x <- which(maksym[1:3358,i] >= 10)
    y <- which(colombos[1:3358,j] >= 10)
    z <- intersect(x,y)
    plot(log10(maksym[z,i]),log10(colombos[z,j]),main="Maksym et al vs COLOMBOS",xlab=bquote("Maksym et al, Sample " ~ .(i-1) ~ ", log10(mRNA abundance)"),ylab=bquote("COLOMBOS, Sample " ~ .(j-1) ~ ", log10(mRNA abundance)"))
    rho.p <- cor(log10(maksym[z,i]),log10(colombos[z,j]))
    xlim <- c(min(log10(maksym[z,i])),max(log10(maksym[z,i])))
    ylim <- c(min(log10(colombos[z,j])),max(log10(colombos[z,j])))
    text(xlim[1]+0.5,ylim[2]-0.25,labels=bquote(rho["p"] ~ "=" ~ .(rho.p)))
    bt.lab.cor.p <- c(bt.lab.cor.p,rho.p)
  }
}
dev.off()

pdf("../Graphs/maksym_li_comp.pdf")
for (i in 2:ncol(maksym))
{
    if ((i-1) %in% maksym.exclude) next
    x <- which(maksym[1:3358,i] >= 10)
    y <- which(li_mrna[1:3358,2] != 0)
    z <- intersect(x,y)
    plot(log10(maksym[z,i]),log10(li_mrna[z,2]),main="Maksym et al vs Li et al",xlab=bquote("Maksym et al, Sample " ~ .(i-1) ~ ", log10(mRNA abundance)"),ylab=bquote("Li et al, log10(mRNA abundance)"))
    rho.p <- cor(log10(maksym[z,i]),log10(li_mrna[z,2]))
    xlim <- c(min(log10(maksym[z,i])),max(log10(maksym[z,i])))
    ylim <- c(min(log10(li_mrna[z,2])),max(log10(li_mrna[z,2])))
    text(xlim[1]+0.5,ylim[2]-0.25,labels=bquote(rho["p"] ~ "=" ~ .(rho.p)))
    bt.lab.cor.p <- c(bt.lab.cor.p,rho.p)
}
dev.off()


pdf("../Graphs/maksym_cho_comp.pdf")
for (i in 2:ncol(maksym))
{
  if ((i-1) %in% maksym.exclude) next
  x <- which(maksym[1:3358,i] >= 10)
  y <- which(cho[1:3358,2] >= 10)
  z <- intersect(x,y)
  plot(log10(maksym[z,i]),log10(cho[z,2]),main="Maksym et al vs Cho et al",xlab=bquote("Maksym et al, Sample " ~ .(i-1) ~ ", log10(mRNA abundance)"),ylab=bquote("Cho et al, log10(mRNA abundance)"))
  rho.p <- cor(log10(maksym[z,i]),log10(cho[z,2]))
  xchom <- c(min(log10(maksym[z,i])),max(log10(maksym[z,i])))
  ychom <- c(min(log10(cho[z,2])),max(log10(cho[z,2])))
  text(xlim[1]+0.5,ylim[2]-0.25,labels=bquote(rho["p"] ~ "=" ~ .(rho.p)))

  bt.lab.cor.p <- c(bt.lab.cor.p,rho.p)
}
dev.off()

pdf("../Graphs/colombos_cho_comp.pdf")
for (i in 2:ncol(colombos))
{
  if ((i-1) %in% colombos.exclude) next
  x <- which(colombos[1:3358,i] >= 10)
  y <- which(cho[1:3358,2] >= 10)
  z <- intersect(x,y)
  plot(log10(colombos[z,i]),log10(cho[z,2]),main="COLOMBOS vs Cho et al",xlab=bquote("COLOMBOS, Sample " ~ .(i-1) ~ ", log10(mRNA abundance)"),ylab=bquote("Cho et al, log10(mRNA abundance)"))
  rho.p <- cor(log10(colombos[z,i]),log10(cho[z,2]))
  xchom <- c(min(log10(colombos[z,i])),max(log10(colombos[z,i])))
  ychom <- c(min(log10(cho[z,2])),max(log10(cho[z,2])))
  text(xlim[1]+0.5,ylim[2]-0.25,labels=bquote(rho["p"] ~ "=" ~ .(rho.p)))
  bt.lab.cor.p <- c(bt.lab.cor.p,rho.p)
}
dev.off()

pdf("../Graphs/colombos_li_comp.pdf")
for (i in 2:ncol(colombos))
{
  if ((i-1) %in% colombos.exclude) next
  x <- which(colombos[1:3358,i] >= 10)
  y <- which(li_mrna[1:3358,2]  != 0)
  z <- intersect(x,y)
  plot(log10(colombos[z,i]),log10(li_mrna[z,2]),main="Colombos et al vs Li et al",xlab=bquote("COLOMBOS, Sample " ~ .(i-1) ~ ", log10(mRNA abundance)"),ylab=bquote("Li et al, log10(mRNA abundance)"))
  rho.p <- cor(log10(colombos[z,i]),log10(li_mrna[z,2]))
  xli_mrnam <- c(min(log10(colombos[z,i])),max(log10(colombos[z,i])))
  yli_mrnam <- c(min(log10(li_mrna[z,2])),max(log10(li_mrna[z,2])))
  text(xlim[1]+0.5,ylim[2]-0.25,labels=bquote(rho["p"] ~ "=" ~ .(rho.p)))
  bt.lab.cor.p <- c(bt.lab.cor.p,rho.p)
}
dev.off()

tmp <- li_mrna[,2] * li_au[,2]
x <- which(tmp[1:3358] != 0)
rho.p <- cor(log10(tmp[x]),phi[x,2])
phi.cor.p <- c(phi.cor.p,rho.p)


tmp <- li_mops_complete[,2]
x <- which(tmp[1:3358] != 0)
rho.p <- cor(log10(tmp[x]),phi[x,2])
phi.cor.p <- c(phi.cor.p,rho.p)

tmp <- li_mops_min[,2]
x <- which(tmp[1:3358] != 0)
rho.p <- cor(log10(tmp[x]),phi[x,2])
phi.cor.p <- c(phi.cor.p,rho.p)

tmp <- li_mops_no_meth[,2]
x <- which(tmp[1:3358] != 0)
rho.p <- cor(log10(tmp[x]),phi[x,2])
phi.cor.p <- c(phi.cor.p,rho.p)

for (i in 2:ncol(maksym))
{

  if (i %in% maksym.exclude) next
  tmp <- maksym[1:3358,i]
  x <- which(tmp >= 10) #exclude genes that appear in the noise
  rho.p <- cor(log10(tmp[x]),phi[x,2])
  phi.cor.p <- c(phi.cor.p,rho.p)
}

for (i in 2:ncol(colombos))
{

  if (i %in% colombos.exclude) next
  tmp <- colombos[1:3358,i]
  x <- which(tmp >= 10) # exclude genes that appear in the noise
  rho.p <- cor(log10(tmp[x]),phi[x,2])
  phi.cor.p <- c(phi.cor.p,rho.p)

}

tmp <- li_mrna[,2]
x <- which(tmp[1:3358] != 0)
rho.p <- cor(log10(tmp[x]),phi[x,2])
phi.cor.p <- c(phi.cor.p,rho.p)

tmp<- cho[,2]
x <- which(tmp[1:3358] >= 10) #exclude genes that appear in the noise
rho.p <- cor(log10(tmp[x]),phi[x,2])
phi.cor.p <- c(phi.cor.p,rho.p)

y <- data.frame(Correlation=c(in.lab.cor.p,bt.lab.cor.p,phi.cor.p),Group=c(rep(1,length(in.lab.cor.p)),rep(2,length(bt.lab.cor.p)),rep(3,length(phi.cor.p))))
y$Group <- factor(y$Group,labels=c("Within Labs","Between Labs","ROC-SEMPPR"))
p <- (ggplot(y,aes(x=Group,y=Correlation))
      +  geom_boxplot(fill="deepskyblue")
      + scale_y_continuous(name=expression(rho~"(Pearson correlation)"),limits=c(0.0,1.0),breaks=seq(0.0,1.0,0.1))
      + scale_x_discrete(name="")
      + theme_bw()
      + theme(axis.title.y=element_text(size = 10,face="bold"),axis.text.x=element_text(size = 10,face = "bold",angle = 45, hjust = 1),axis.text.y=element_text(size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title=element_text(size=12,hjust=0.5))
      + ggtitle(label="Comparison of Gene Expression\nMeasurements (Endogenous only)"))
ggsave("../Graphs/boxplots_test.pdf",device="pdf",dpi = 300)