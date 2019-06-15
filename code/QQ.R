#!/usr/bin/env Rscript
library("argparse")
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--pvalfile", action="store")
parser$add_argument("--plotfile", action="store")
parser$add_argument("--maxp", default=100, action="store", type="integer")
parser$add_argument("--recessive", default=FALSE, action="store_true")
parser$add_argument("--removenull", default=FALSE, action="store_true")
args <- parser$parse_args()

dat<-read.delim(args$pvalfile, header=T, stringsAsFactors = F, sep="\t")
plotname<-args$plotfile

if(args$recessive==TRUE){
  dat$log10p_obs<-(-log10(dat$P_REC))
}else{
  dat$log10p_obs<-(-log10(dat$P_DOM))
}

if(args$removenull==TRUE){
  dat<-subset(dat, log10p_obs>0)
}

dat<-dat[order(dat$log10p_obs),]
dat$log10p_exp<-sort(-log10(runif(nrow(dat))))

#Calculate slope of line 
p0<-0
u0=max(dat[dat$log10p_obs==0,]$log10p_exp)
tmp<-subset(dat, log10p_obs!=0)
p95<-tmp[round(nrow(tmp[tmp$log10p_obs!=0,]))*0.95,]$log10p_obs
u95<-tmp[round(nrow(tmp[tmp$log10p_obs!=0,]))*0.95,]$log10p_exp
slope<-(p95-p0)/(u95-u0)
yint<-p95-slope*u95

maxp<-ceiling(max(dat$log10p_obs))
if(maxp>args$maxp){maxp<-args$maxp}

#Print Lambda
print(paste("lambda=", slope, sep=""))

#Generate the actual plot
options(bitmapType='cairo')
png(plotname, width=500, height=500)
par(mar=c(4,5,2,2))
plot(x=dat$log10p_exp, y=dat$log10p_obs, ylim=c(0, maxp), xlim=c(0,maxp), cex=0.5, pch=19,cex.lab=1.5, cex.axis=1.3,xlab=expression(Expected ~ -log[10](p)), ylab=expression(Observed ~ -log[10](p)), main="")
abline(yint, slope, col="blue", lty=2)
abline(0,1)
dev.off()
