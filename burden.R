#!/usr/bin/env Rscript
library("argparse")
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--casefile", action="store")
parser$add_argument("--casesize", action="store", type="integer")
parser$add_argument("--controlfile", action="store")
parser$add_argument("--controlsize", action="store", type="integer")
parser$add_argument("--outfile", action="store")

args <- parser$parse_args()
#args$casesize+args$controlsize

case.dat<-read.delim(args$casefile, header=T, stringsAsFactors=F, sep="\t")
names(case.dat)[1]<-"GENE"
control.dat<-read.delim(args$controlfile, header=T, stringsAsFactors=F, sep="\t")
names(control.dat)[1]<-"GENE"

dat<-merge(case.dat, control.dat, by="GENE", all.x=T, all.y=T)
dat[is.na(dat)]<-0

dat$P_DOM<-0
dat$P_REC<-0

for(i in 1:nrow(dat)){
  mat<-cbind(c(dat[i,]$CASE_COUNT_ALL, (args$casesize-dat[i,]$CASE_COUNT_ALL)), c(dat[i,]$CONTROL_AC_ALL, (args$controlsize-dat[i,]$CONTROL_AC_ALL)))
  dat[i,]$P_DOM<-fisher.test(mat, alternative="greater")$p.value
}

write.table(x=dat,file=args$outfile, sep="\t", quote=F, row.names=F, col.names=T)
