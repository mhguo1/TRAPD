#!/usr/bin/python
import optparse
import sys
import gzip 

#Parse options
parser = optparse.OptionParser()

#Required Options
parser.add_option("-s", "--snpfile", action="store",dest="snpfilename")
parser.add_option("-v", "--vcffile", action="store",dest="vcffilename")
parser.add_option("-o", "--outfile", action="store",dest="outfilename", default="out.txt")

parser.add_option("--snpformat", action="store",dest="snpformat", default="VCFID")
parser.add_option("--pop", action="append",dest="pop", default="ALL")

#Optional Filters
parser.add_option("--pass", action="store_true", dest="passfilter")
parser.add_option("--maxAF", action="store",dest="maxAF", default=1)
parser.add_option("--maxAC", action="store",dest="maxAC", default=99999)

options, args = parser.parse_args()
if not options.snpfilename:   # if filename is not given
    parser.error('A file containing a list of SNPs is needed')

if not options.vcffilename:   # if filename is not given
    parser.error('A vcf file is needed')


def makesnplist(snpfile):
	#Makes a list of SNPs present in the snpfile
	snplist=[]
	#Read in snpfile
	snp_file=open(snpfile, "r")
	
	for line_snp1 in snp_file:
		line_snp=line_snp1.rstrip('\n').split('\t')

		#Find column corresponding to desired snps
		if line_snp[0]!="GENE":
			snplist=snplist+line_snp[1].split(",")
	return set(snplist)
	snp_file.close()

def extractcounts(pops, vcfline):
	ac_out=0
	ac_hom_out=0
	af_out=0
	if (options.pop is None) or ("ALL" in options.pop):
		ac_out=(";"+vcfline).split((";AC="))[1].split(";")[0]
		ac_hom_out=(";"+vcfline).split((";Hom="))[1].split(";")[0]
	else:
		for p in range(0, len(pops), 1):
			temp_pop=pops[p]
			ac_out=ac_out+int((";"+vcfline).split((";AC"+pop+"="))[1].split(";")[0])
			ac_hom_out=ac_hom_out+int((";"+vcfline).split((";Hom_"+pop+"="))[1].split(";")[0])
	return [ac_out, ac_hom_out]


def sumcount(snplist, snptable):
	ac_sum=0
	ac_hom_sum=0
	for s in range(0, len(snplist), 1):
		if snplist[s] in snptable:
			tempsnp=snplist[s]
			ac_sum=ac_sum+int(snptable[tempsnp][1])
			ac_hom_sum=ac_hom_sum+int(snptable[tempsnp][2])
	return [ac_sum, ac_hom_sum]

#Make list of all SNPs across all genes present in snpfile
allsnplist=makesnplist(options.snpfilename)

#Make a hashtable with keys as each SNP, and stores a list of indices of carriers for that SNP
count_table={} 

#Open vcf file
if str(options.filename)[-3:]==".gz":
	vcffile=gzip.open(options.vcffilename, "rb")
else:
	vcffile=open(options.vcffilename, "r")
for line_vcf1 in vcffile:
	line_vcf=line_vcf1.rstrip().split('\t')
	if line_vcf[0][0]!="#":
		if not (options.passfilter and line_vcf[6]!="PASS"):
			if options.snpformat=="VCF":
				snpid=str(line_vcf[2])
			else: 
				snpid=str(line_vcf[0].lstrip("chr"))+":"+str(line_vcf[1])+":"+str(line_vcf[3])+":"+str(line_vcf[4])
			if snpid in listofsnps:
				counts=extractcounts(options.pop, line_vcf[7])
				count_table[snpid]=[snpid, counts[0], counts[1]]
vcffile.close()


#Test by writing output
outfile=open(options.outfilename, "w")
snpfile=open(options.snpfilename, "r")
for line_s1 in snpfile:
	line_s=line_s1.rstrip('\n').split('\t')
	if line_s[0]!="GENE":
		genesnplist=list(set(line_s[1].split(',')))
		sumcounts=sumcount(genesnplist, count_table)
		outfile.write(line_s[0]+"\t"+str(sumcounts[0])+"\t"+str(sumcounts[1])+'\n')
outfile.close()
snpfile.close()


