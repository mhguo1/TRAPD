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

case.dat<-read.delim(args$casefile, header=T, stringsAsFactors=F, sep="\t")
names(case.dat)[1]<-"GENE"
control.dat<-read.delim(args$controlfile, header=T, stringsAsFactors=F, sep="\t")
names(control.dat)[1]<-"GENE"

dat<-merge(case.dat, control.dat, by="GENE", all.x=T, all.y=T)
dat[is.na(dat)]<-0

dat$P_DOM<-0
dat$P_REC<-0

for(i in 1:nrow(dat)){
  
  #Dominant model
  case_count<-dat[i,]$CASE_COUNT_HET+dat[i,]$CASE_COUNT_CH+dat[i,]$CASE_COUNT_HOM
  control_count<-dat[i,]$CONTROL_COUNT_HET+dat[i,]$CONTROL_COUNT_HOM
  
  if(case_count>args$casesize){
    case_count<-args$casesize
  }else if(case_count<0){
    case_count<-0
   }
  if(control_count>args$controlsize){
    control_count<-args$controlsize
  }else if(control_count<0){
    control_count<-0
   }
  
  mat<-cbind(c(case_count, (args$casesize-case_count)), c(control_count, (args$controlsize-control_count)))
  dat[i,]$P_DOM<-fisher.test(mat, alternative="greater")$p.value
  
  
  #Recessive model
  case_count_rec<-dat[i,]$CASE_COUNT_CH+dat[i,]$CASE_COUNT_HOM
  control_count_rec<-dat[i,]$CONTROL_COUNT_HOM+(args$controlsize)*((dat[i,]$CONTROL_COUNT_HET)/(args$controlsize))^2
  
  if(control_count_rec<0){ control_count_rec<-0}
  if(case_count_rec>args$casesize){case_count_rec<-args$casesize}
  if(control_count_rec>args$controlsize){control_count_rec<-args$controlsize}
  control_count_rec<-round(control_count_rec, digits=0)
  
  mat_rec<-cbind(c(case_count_rec, (args$casesize-case_count_rec)), c(control_count_rec, (args$controlsize-control_count_rec)))
  dat[i,]$P_REC<-fisher.test(mat_rec, alternative="greater")$p.value
}

write.table(x=dat,file=args$outfile, sep="\t", quote=F, row.names=F, col.names=T)
