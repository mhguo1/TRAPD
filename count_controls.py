#!/usr/bin/python
import optparse
import sys
import gzip 

#Parse options
parser = optparse.OptionParser()
parser.add_option("--snpfile", action="store",dest="snpfilename")
parser.add_option("--snpformat", action="store",dest="snpformat", default="VCFID")
parser.add_option("--snpcolname", action="store",dest="snpcolname", default="NA")
parser.add_option("--controlvcf", action="store",dest="vcffilename")
parser.add_option("--pop", action="append",dest="pop", default="ALL")
parser.add_option("--outfile", action="store",dest="outfilename")
parser.add_option("--recessive", action="store_true",dest="recessive")
parser.add_option("--maxAF", action="store",dest="maxAF", default=1)
parser.add_option("--maxAC", action="store",dest="maxAC", default=99999)

options, args = parser.parse_args()
if not options.snpfilename:   # if filename is not given
    parser.error('A file containing a list of SNPs is needed')

if not options.vcffilename:   # if filename is not given
    parser.error('A vcf file is needed')


def makesnplist(snpfile, snpcolname):
	
	snplist=[]

	#Read in snpfile
	snp_file=open(snpfile, "r")
	
	for line_snp1 in snp_file:
		line_snp=line_snp1.rstrip('\n').split('\t')

		#Find column corresponding to desired snps
		if line_snp[0]=="GENE":
			if snpcolname!="NA":
				snpcol=line_snp.index(snpcolname)
			else:
				snpcol=1
		elif len(line_snp[snpcol])>1:
			#Add SNPs to list
			snplist=snplist+line_snp[snpcol].split(",")

	return snplist
	snp_file.close()


def extractcounts(pops, vcfline):
	ac_out=0
	af_out=0
	if (options.pop is None) or ("ALL" in options.pop):
		
		ac_out=(";"+vcfline).split((";AC="))[1].split(";")[0]
		af_out=(";"+vcfline).split((";AF="))[1].split(";")[0]
	else:
		for p in range(0, len(pops), 1):
			temp_pop=pops[p]
			ac_out=ac_out+int((";"+vcfline).split((";AC"+pop+"="))[1].split(";")[0])
			af_out=af_out+float((";"+vcfline).split((";AF"+pop+"="))[1].split(";")[0])
	return [ac_out, af_out]


def sumcount(snplist, snptable):
	ac_sum=0
	af_sum=0
	for s in range(0, len(snplist), 1):
		if snplist[s] in snptable:
			tempsnp=snplist[s]
			ac_sum=ac_sum+int(snptable[tempsnp][1])
			af_sum=af_sum+float(snptable[tempsnp][2])
	return [ac_sum, af_sum]




listofsnps=makesnplist(options.snpfilename, options.snpcolname)
tableout={}
#Open vcf file
if str(options.filename)[-3:]==".gz":
	vcffile=gzip.open(options.vcffilename, "rb")
else:
	vcffile=open(options.vcffilename, "r")

for line_vcf1 in vcffile:
	line_vcf=line_vcf1.rstrip().split('\t')
	if line_vcf[0][0]!="#":
		if options.snpformat=="VCF":
			snpid=str(line_vcf[2])
		else: 
			snpid=str(line_vcf[0].lstrip("chr"))+":"+str(line_vcf[1])+":"+str(line_vcf[3])+":"+str(line_vcf[4])
		if snpid in listofsnps:
			line_out=extractcounts(options.pop, line_vcf[7])
			tableout[snpid]=[snpid, line_out[0], line_out[1]]

vcffile.close()


#Test by writing output
outfile=open(options.outfilename, "w")
snpfile=open(options.snpfilename, "r")
for line_s1 in snpfile:
	line_s=line_s1.rstrip('\n').split('\t')
	if line_s[0]!="GENE":
		#snplist=list(set(line_s[1].split(',')))
		snplist=line_s[1].split(',')
		sumcounts=sumcount(snplist, tableout)
		outfile.write(line_s[0]+"\t"+str(sumcounts[0])+"\t"+str(sumcounts[1])+'\n')
outfile.close()
snpfile.close()


#python count_control.py --controlvcf gnomad.test.vcf.gz --snpfile test.out.txt --out test.control.out.txt
