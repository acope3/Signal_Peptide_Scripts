fr <- function(x,b,a,x.obs,y.obs,sd.1,sd.2)
{
  y <- b*x + a
  dnorm(x.obs,mean = x,sd = sd.1) * dnorm(y.obs,mean=y,sd=sd.2)
}

test <- function(b,a)
{
  likelihood <- 0
  for (i in 1:length(x.obs))
  {
    x <- integrate(fr,lower=-Inf,upper=Inf,b=b,a=a,x.obs=x.obs[i],y.obs=y.obs[i],sd.1=sd.1[i],sd.2=sd.2[i],stop.on.error=F)
    likelihood <- likelihood + log(x$value)
  }
  -1*likelihood
}

plotResults<-function(data,b1,b0,categories=c("X","Y"),ci = F,bounds=NULL,file="title",title="Regression",mark=NULL)
{
  
  if (!is.null(bounds))
  {
    conf.int.1 <- bounds[1]
    conf.int.2 <- bounds[2]
    l <- data.frame(s=c(b1,conf.int.2,conf.int.1,1.0),ic=c(b0,0.0,0.0,0.0),Line=c("Model II Regression","97.5% CI","2.5% CI","y=x"),stringsAsFactors = F)
    l$Line <- factor(l$Line,levels=c("Model II Regression","97.5% CI","2.5% CI","y=x"))
    levels(l$Line) <- c("Model II Regression","95% CI","95% CI","y=x")
    legend.colors <- c("black","grey","black")
    lines.reg <- c("solid","dashed","dashed")

  } else{
    l <- data.frame(s=c(b1,1.0),ic=c(b0,0.0),Line=c("Model II Regression","y = x"))
    legend.colors <- c("black","black")
    lines.reg <- c("solid","dashed")
  }
  p <- ggplot(data,aes(Posterior,Posterior.2))
  p <-(p + geom_point(colour="black",size=4)
       + labs(x=bquote(.(categories[1])~"("*Delta*eta*")"),y=bquote(.(categories[2])~"("*Delta*eta*")")) 
       + geom_abline(data=l,mapping=aes(slope=s,intercept=ic,linetype=Line,color=Line))
       + scale_color_manual(values=legend.colors)
       + scale_linetype_manual(values=lines.reg)
       + ggtitle(label=title))
  p<- p + geom_text(aes(label=ifelse(data$Codon %in% mark,as.character(data$Codon),'')),nudge_x =-0.2,nudge_y=0.1)
  
  rho.p <- round(cor(data[,"Posterior"],data[,"Posterior.2"]),2)
  if (ci)
  {
    p <- (p + geom_errorbar(mapping=aes(ymin=X0.025.2,ymax=X0.975.2),alpha=0.2,color="black") 
          + geom_errorbarh(mapping=aes(xmin=X0.025.,xmax=X0.975.),alpha=0.2,color="black"))
    range.xy <- range(c(data[,4:5],data[,7:8]),na.rm=T)
    xlim <- range.xy
    ylim <- range.xy
  } else{
    range.xy <- range(c(data[,3],data[,6]),na.rm=T)
    xlim <- range.xy
    ylim <- range.xy
  }
  p <- p + scale_x_continuous(limits = range.xy) + scale_y_continuous(limits = range.xy)
  width <- xlim[2] - xlim[1]
  height <- ylim[2] - ylim[1]
  cor.exp <- bquote(rho ~ " = " ~ .(rho.p))
  b1 <- round(b1,3)
  b0 <- round(b0,3)
  if (b0 != 0)
  {
    eq.exp <- bquote("y ="~ .(format(b1,nsmall=3))*"x +"~ .(format(b0,nsmall=3)))
  } else{
    eq.exp <- bquote("y ="~ .(format(b1,nsmall=3))*"x")
  }
  p <- p + annotate("text",x=xlim[2] - width * 0.2,y=ylim[1] + height * 0.10,label=deparse(cor.exp),parse=T,size=5) 
  p <- p + annotate("text",x=xlim[2]-width*0.2,y=ylim[1]+height*0.20,label=deparse(eq.exp),parse=T,size=5)
  p <- (p + theme_bw()
        + theme(axis.title=element_text(size = 12,face="bold"),axis.text=element_text(size=12))
        + theme(axis.line = element_line(colour = "black"))
        + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        + theme(legend.position = c(0.25,0.8),legend.title=element_blank(),legend.text=element_text(size=12),plot.title = element_text(hjust = 0.5,size=12)))
  
  
  ggsave(filename = file,device="pdf",width = 7,height = 7)
  return(p)
}

normalize <- function(dEta)
{
  aa <- AnaCoDa::aminoAcids()
  for (a in aa)
  {
    if (a=="W" || a=="M"||a=="X") next
    row.ind <- which(dEta[,1] == a)
    rows <- dEta[row.ind,]
    mu <- mean(rows[,3])
    rows[,3:5] <- rows[,3:5] - mu
    dEta[row.ind,] <- rows
  }
  return(dEta)
}

getConfInt<-function(dEta)
{
  aa <- AnaCoDa::aminoAcids()
  for (a in aa)
  {
    if (a=="W" || a=="M"||a=="X") next
    row.ind <- which(dEta[,1] == a)
    rows <- dEta[row.ind,]
    index <- which.max(rows[,5] - rows[,4])
    #index <- sample(1:length(row.ind),size = 1)
    ref <- which(rows[,3] == 0)
    rows[ref,4] <- 0 - abs(rows[index,4] - rows[index,3]) 
    rows[ref,5] <- 0 + abs(rows[index,5] - rows[index,3])
    dEta[row.ind,] <- rows
  }
  return(dEta)
}

rescaleCSPEstimates<- function(param.1)
{
  param.new <- data.frame()
  aa <- unique(param.1[,1])
  for (a in aa)
  {
    codons <- AAToCodon(a)
    tmp <- param.1[codons,]
    na.row <- which(is.na(tmp[,1]))
    rownames(tmp)[na.row] <- codons[na.row]
    tmp[na.row,1:2] <-c(a,codons[na.row])
    tmp[na.row,3:5] <-c(0,0,0)
    param.new <- rbind(param.new,tmp)
  }
  return(param.new)
}

fixParameter <- function(parameter.file.to.fix)
{
  load(parameter.file.to.fix)
  
  aminoAcids()
  paramBase$grouplist <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R", "S", "T", "V", "Y", "Z")
  
  save(list = c("paramBase", "currentMutation", "currentSelection",
                "proposedMutation", "proposedSelection", "model",  
                "mutationPrior", "mutationTrace", "selectionTrace", 
                "synthesisOffsetAcceptRatTrace", "synthesisOffsetTrace", 
                "observedSynthesisNoiseTrace", "withPhi"),
       file=parameter.file.to.fix)
}

library(AnaCoDa)
library(ggplot2)
library(bbmle)
directory <- "../Results/mp_sp_pseudo_emp_imputed/chain/run_4/"
parameter <- loadParameterObject(paste0(directory,"R_objects/parameter.Rda"))
names.aa <- aminoAcids()
exclude <- c()
sd.1 <- vector(mode = "numeric",length=59 - length(exclude))
sd.2 <- vector(mode="numeric",length=59 - length(exclude))
index <- 1
for(aa in names.aa)
{
  if(aa == "M" || aa == "W" || aa == "X") next
  codons <- AAToCodon(aa,T)
  for(j in 1:length(codons))
  {
    if (codons[j] %in% exclude) next
    sd.1[index] <- parameter$getCodonSpecificVariance(codons[j],4, 20000, 1, TRUE, TRUE)
    sd.2[index] <- parameter$getCodonSpecificVariance(codons[j],2, 20000, 1, TRUE, TRUE)
    index <- index + 1
  }
  sd.1[index] <- sum(sd.1[(index-length(codons)):(index-1)])/(length(codons)+1)
  sd.2[index] <- sum(sd.2[(index-length(codons)):(index-1)])/(length(codons)+1)
  sd.1[(index-length(codons)):(index-1)] <- sqrt(sd.1[index])
  sd.2[(index-length(codons)):(index-1)] <- sqrt(sd.2[index])
  sd.1[index] <- sqrt(sd.1[index])
  sd.2[index] <- sqrt(sd.2[index])
  index <- index +1
}
rm(parameter)
sel.1 <- read.table(paste0(directory,"Parameter_est/pseudo_sp_w_emp_final_Selection"),sep=",",header=TRUE,stringsAsFactors=F)
sel.2 <- read.table(paste0(directory,"Parameter_est/sp_main_Selection"),sep=",",header=TRUE,stringsAsFactors=F)
rownames(sel.1) <- sel.1[,2]
rownames(sel.2) <- sel.2[,2]
#sel.1 <- rescaleCSPEstimates(sel.1)
#sel.2 <- rescaleCSPEstimates(sel.2)
sel.1 <- getConfInt(sel.1)
sel.2 <- getConfInt(sel.2)
sel.1 <- normalize(sel.1)
sel.2 <- normalize(sel.2)
#sel.1 <- sel.1[-c(which(sel.1[,2] %in% exclude)),]
#sel.2 <- sel.2[-c(which(sel.2[,2] %in% exclude)),]
reg <- lm(sel.2[,3] ~ sel.1[,3]+0)
#line.values <- mle2(test,start=list(b=reg$coefficients[[2]],a=reg$coefficients[[1]]),data=list(x.obs=sel.1[,3],y.obs=sel.2[,3],sd.1=sd.1,sd.2=sd.2,lower.x=sel.1[,4],upper.x=sel.1[,5]))
line.values <- mle2(test,start=list(b=reg$coefficients[[1]]),fixed=list(a=0),data=list(x.obs=sel.1[,3],y.obs=sel.2[,3],sd.1=sd.1,sd.2=sd.2,lower.x=sel.1[,4],upper.x=sel.1[,5]))
prof <- profile(line.values)
#line.values <- read.table("../../../For_paper/psp_sp_51_no_int.txt",sep="",header=F)[,1]
df <- plyr::join(sel.1,sel.2,by=c("AA","Codon"))
colnames(df)[6] <- "Posterior.2"
colnames(df)[7] <- "X0.025.2"
colnames(df)[8] <- "X0.975.2"
ci <- confint(prof)
p<-plotResults(data=df,b1=line.values@coef[[1]],b0=0,ci=T,bounds=ci,file="../../../For_paper/check_w_ci.pdf",categories=c("Mature Peptides","Signal Peptides"),title = "Selection on Codons:\nMature vs. Signal Peptides")
plotResults(data=df,b1=line.values@coef[[1]],b0=0,ci=F,bounds=ci,file="../../../For_paper/check_wo_ci.pdf",categories=c("Mature Peptides","Signal Peptides"),title = "Selection on Codons:\nMature vs. Signal Peptides")
reg_result <- c(line.values@coef[[1]],0,ci[1],ci[2])
write(reg_result,"~/For_paper/check.txt",ncolumns = 1)

# slopes<-c()
# lower.bounds <- c()
# upper.bounds <- c()
# corr <- c()
# for(i in 2:35)
# {
#   cat(i,"\n")
#   directory <-paste0("../Results/Segmented_Genome_main/wphi/Segment_1_",i,"/chain",(i-1),"/run_2/")
#   fixParameter(paste0(directory,"R_objects/parameter.Rda"))
#   parameter <- loadParameterObject(paste0(directory,"R_objects/parameter.Rda"))
#   names.aa <- aminoAcids()
#   exclude <- c()
#   sd.1 <- vector(mode = "numeric",length=59 - length(exclude))
#   sd.2 <- vector(mode="numeric",length=59 - length(exclude))
#   index <- 1
#   for(aa in names.aa)
#   {
#     if(aa == "M" || aa == "W" || aa == "X") next
#     codons <- AAToCodon(aa,T)
#     for(j in 1:length(codons))
#     {
#       if (codons[j] %in% exclude) next
#       sd.1[index] <- parameter$getCodonSpecificVariance(codons[j],1, 5000, 1, TRUE, TRUE)
#       sd.2[index] <- parameter$getCodonSpecificVariance(codons[j],2, 5000, 1, TRUE, TRUE)
#       index <- index + 1
#     }
#     sd.1[index] <- sum(sd.1[(index-length(codons)):(index-1)])/(length(codons)+1)
#     sd.2[index] <- sum(sd.2[(index-length(codons)):(index-1)])/(length(codons)+1)
#     sd.1[(index-length(codons)):(index-1)] <- sqrt(sd.1[index])
#     sd.2[(index-length(codons)):(index-1)] <- sqrt(sd.2[index])
#     sd.1[index] <- sqrt(sd.1[index])
#     sd.2[index] <- sqrt(sd.2[index])
#     index <- index +1
#   }
#   rm(parameter)
#   sel.1 <- read.table(paste0(directory,"Parameter_est/selection_segment_1_main.csv"),sep=",",header=TRUE,stringsAsFactors = F)
#   sel.2 <- read.table(paste0(directory,"Parameter_est/selection_segment_",i,"_main.csv"),sep=",",header=TRUE,stringsAsFactors = F)
#   rownames(sel.1) <- sel.1[,2]
#   rownames(sel.2) <- sel.2[,2]
#   sel.1 <- rescaleCSPEstimates(sel.1)
#   sel.2 <- rescaleCSPEstimates(sel.2)
#   sel.1 <- getConfInt(sel.1)
#   sel.2 <- getConfInt(sel.2)
#   sel.1 <- normalize(sel.1)
#   sel.2 <- normalize(sel.2)
#   #sel.1 <- sel.1[-c(which(sel.1[,2] %in% exclude)),]
#   #sel.2 <- sel.2[-c(which(sel.2[,2] %in% exclude)),]
#   reg <- lm(sel.2[,3] ~ sel.1[,3]+0)
#   line.values <- mle2(test,start=list(b=reg$coefficients[[1]]),fixed=list(a=0),data=list(x.obs=sel.1[,3],y.obs=sel.2[,3],sd.1=sd.1,sd.2=sd.2,lower.x=sel.1[,4],upper.x=sel.1[,5]),method="SANN")
#   prof <- profile(line.values)
#   ci <- confint(prof)
#   slopes<-c(slopes,line.values@coef[[1]])
#   lower.bounds <- c(lower.bounds,ci[1])
#   upper.bounds <- c(upper.bounds,ci[2])
#   corr <- c(corr,cor(sel.2[,3],sel.1[,3]))
# }
# x<-c(2:35)
# pdf(paste0("CUB_variation_Ecoli_genome_ROC_rev.pdf"))
# par(mar=c(5,5,5,5))
# plot(x,slopes,xlab="Segment",ylab="Regression Slope",xlim = c(2,35),main=paste0("CUB Variation"),ylim=c(min(lower.bounds)-0.05,max(upper.bounds)+0.05))
# arrows(x,lower.bounds, x, upper.bounds, length=0.05, angle=90, code=3)
# par(new=TRUE)
# plot(x,corr,pch=17, axes=F,xlab=NA,ylab=NA,cex=1.2,ylim=c(min(corr)-0.1,max(corr)+0.1))
# axis(side=4)
# mtext(side=4,line=3,"Pearson Correlation")
# legend("topright",inset=0.05,legend = c("Slope","Pearson Correlation"),pch=c(1,17))
# dev.off()
# df <- data.frame(slope=slopes,lower.bound=lower.bounds,upper.bound=upper.bounds,correlation=corr)
# write.table(df,file = "Ecoli_K12_genome_ROC_rev.csv",sep=",",row.names = F,quote=F)

