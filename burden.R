#!/usr/bin/env Rscript
library("argparse")
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--casefile", action="store")
parser$add_argument("--casesize", action="store", type="integer")
parser$add_argument("--controlfile", action="store")
parser$add_argument("--controlsize", action="store", type="integer")

args <- parser$parse_args()
args$casesize+args$controlsize

case.dat<-read.delim(args$casefile, header=T, stringsAsFactors=F, sep="\t")
names(case.dat)[1]<-"GENE"
control.dat<-read.delim(args$controlfile, header=T, stringsAsFactors=F, sep="\t")
names(control.dat)[1]<-"GENE"

dat<-merge(case.dat, control.dat, by="GENE", all.x=T, all.y=T)
dat[is.na(dat)]<-0



##Rscript test.R --casesize 2 --controlsize 3
