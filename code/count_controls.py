#!/usr/bin/python
import optparse
import operator
import re
import sys
import gzip 
import pybedtools
from pybedtools import BedTool

#Parse options
parser = optparse.OptionParser()

#Required Options
parser.add_option("-s", "--snpfile", action="store",dest="snpfilename")
parser.add_option("-v", "--vcffile", action="store",dest="vcffilename")

parser.add_option("-o", "--outfile", action="store",dest="outfilename", default="control_counts.txt")

parser.add_option("--snpformat", action="store",dest="snpformat", default="VCFID")
parser.add_option("--pop", action="store",dest="pop", default="ALL")
parser.add_option("-d", "--database", action="store", dest="database", default="generic")

#Optional Filters
parser.add_option("--pass", action="store_true", dest="passfilter")
parser.add_option("--maxAF", action="store",dest="maxAF", default=1)
parser.add_option("--maxAC", action="store",dest="maxAC", default=99999)
parser.add_option("--minAN", action="store",dest="minAN", default=0)
parser.add_option("--popmaxAF", action="store",dest="popmaxAF", default=1)
parser.add_option("--homcol", action="store",dest="homcol")
parser.add_option("--bedfile", action="store", dest="bedfilename")

options, args = parser.parse_args()
if not options.snpfilename:   # if filename is not given
    parser.error('A file containing a list of SNPs is needed')

if not options.vcffilename:   # if filename is not given
    parser.error('A vcf file is needed')

if options.database not in ["generic", "gnomad", "exac"]:
	parser.error('Database must be generic, gnomad, or exac')
	
#Parse populations
if options.database!="generic" and options.pop is not None:
	pops=str(options.pop).split(',')
else:
	pops=["ALL"]

#Check to make sure all populations listed are correct
if options.database=="gnomad":
	pop_list=["AFR", "AMR", "ALL", "ASJ", "EAS", "FIN", "NFE", "SAS"]
	if not all(p in pop_list for p in pops):
		   parser.error('Please check the populations listed')
if options.database=="exac":
	pop_list=["AFR", "AMR", "ALL", "EAS", "FIN", "NFE", "SAS"]
	if not all(p in pop_list for p in pops):
		   parser.error('Please check the populations listed')
if options.database=="generic" and "ALL" not in pops:
	parser.error('Please check the populations listed for your control database')

#Find out chromosome format
vcffile=gzip.open(options.vcffilename, "rb")
chrformat="number"
for line_vcf1 in vcffile:
	line_vcf=line_vcf1.split("\t")
	if "##contig" in line_vcf[0]:
		if "ID=chr" in line_vcf[0]:
			chrformat="chr"
vcffile.close()
			
#Check generic database to make sure it has AC, AN and AF
if options.database=="generic":
	vcffile=gzip.open(options.vcffilename, "rb")
        ac_found=0
	an_found=0
        for line_vcf1 in vcffile:
                if line_vcf1[0]=="#" and ("ID=AC," in line_vcf1):
                        ac_found=1
		elif line_vcf1[0]=="#" and ("ID=AN," in line_vcf1):
                        an_found=1
		elif line_vcf1[0:1]=="#C":
                        break
        if ac_found==0 or an_found==0:
                sys.stdout.write("AC and AN not found in vcf file\n")
                sys.exit()
        vcffile.close()
	
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

def extractcounts(pops, vcfline, max_ac, max_af, popmax_af,min_an):
	ac_out=float((";"+vcfline).split((";AC="))[1].split(";")[0])
	an=float((";"+vcfline).split((";AN="))[1].split(";")[0])
	if ac_out==0:
		af_out=0
	else:
		if options.database=="generic":
			if ";AF=" in vcfline:
				af_out=float((";"+vcfline).split((";AF="))[1].split(";")[0])
			else:
				af_out=float(ac_out)/float(an_out)
		else:	
			af_out=float((";"+vcfline).split((";AF="))[1].split(";")[0])
			
	if popmax_af<1:
		af_popmax_out=get_popmax(vcfline)
	else:
                af_popmax_out=1
		
	        af_popmax_out=float(af_popmax_out)

	if (ac_out>float(max_ac)) or (af_out>float(max_af)) or (an<float(min_an)) or (af_popmax_out>float(popmax_af)):
		ac_out=0
		hom_out=0
	else:
		if options.database=="gnomad":
			if "ALL" in pops:
				hom_out=(";"+vcfline).split((";Hom="))[1].split(";")[0]
			elif "ALL" not in pops:
				ac_out=0
               			hom_out=0
				for p in range(0, len(pops), 1):
					temp_pop=pops[p]
					ac_out=int(ac_out)+int((";"+vcfline).split((";AC_"+temp_pop+"="))[1].split(";")[0])
					hom_out=int(hom_out)+int((";"+vcfline).split((";Hom_"+temp_pop+"="))[1].split(";")[0])
		elif options.database=="exac":
			if "ALL" in pops:
				hom_out=(";"+vcfline).split((";AC_Hom="))[1].split(";")[0]
			elif "ALL" not in pops:
				ac_out=0
               			hom_out=0
				for p in range(0, len(pops), 1):
					temp_pop=pops[p]
					ac_out=int(ac_out)+int((";"+vcfline).split((";AC_"+temp_pop+"="))[1].split(";")[0])
					hom_out=int(hom_out)+int((";"+vcfline).split((";Hom_"+temp_pop+"="))[1].split(";")[0])

		elif options.database=="generic":
			if options.homcol is not None:
				hom_out=(";"+vcfline).split((options.homcol))[1].split(";")[0]
			else:
				hom_out=0	
	return [ac_out, hom_out]


def get_popmax(vcfline):
	if options.database in ["gnomad", "generic"]:
		if ";AF_POPMAX=" in (";"+vcfline):
			af_popmax_out=(";"+vcfline).split((";AF_POPMAX="))[1].split(";")[0]
        		if af_popmax_out==".":
              	 		af_popmax_out=0
		else:
			af_popmax_out=1
	if options.database=="exac":
		if (";AC_POPMAX=" in (";"+vcfline)) and (";AN_POPMAX" in (";"+vcfline)):
			ac_popmax=(";"+vcfline).split((";AC_POPMAX="))[1].split(";")[0]
			an_popmax=(";"+vcfline).split((";AN_POPMAX="))[1].split(";")[0]
			if str(ac_popmax)=="NA" or str(an_popmax)=="NA":
				af_popmax_out=0
			else:
				af_popmax_out=float(ac_popmax)/float(an_popmax)
		else:
			af_popmax_out=1

def sumcount(genesnps, snptable):
	ac_sum=0
	hom_sum=0
	for s in range(0, len(genesnps), 1):
		if genesnps[s] in snptable:
			tempsnp=genesnps[s]
			ac_sum=ac_sum+int(snptable[tempsnp][1])
			hom_sum=hom_sum+int(snptable[tempsnp][2])
	return [ac_sum, hom_sum]

#Make list of all SNPs across all genes present in snpfile
allsnplist=makesnplist(options.snpfilename)

#Make a hashtable with keys as each SNP, and stores a list of indices of carriers for that SNP
count_table={} 

#Open vcf file
vcffile=BedTool(options.vcffilename)
if options.bedfilename is not None:
        bed=BedTool(options.bedfilename)
        vcffile_temp=vcffile.intersect(bed)
else:
        if chrformat=="chr":
                dummy_bed=BedTool('chr1000 100000000 100000001', from_string=True)
        else:
                dummy_bed=BedTool('1000 100000000 100000001', from_string=True)
        vcffile_temp=vcffile.subtract(dummy_bed)

for line_vcf1 in open(vcffile_temp.fn):
	line_vcf=line_vcf1.rstrip().split('\t')
	if line_vcf[0][0]!="#":
		if not (options.passfilter and line_vcf[6]!="PASS"):
			if options.snpformat=="VCFID":
				snpid=str(line_vcf[2])
			else: 
				snpid=str(line_vcf[0].lstrip("chr"))+":"+str(line_vcf[1])+":"+str(line_vcf[3])+":"+str(line_vcf[4])
			if snpid in allsnplist:
				counts=extractcounts(pops, line_vcf[7], options.maxAC, options.maxAF, options.popmaxAF,options.minAN)
				count_table[snpid]=[snpid, counts[0], counts[1]]
pybedtools.cleanup() 

#Write output
outfile=open(options.outfilename, "w")
outfile.write("#GENE\tCONTROL_COUNT_HET\tCONTROL_COUNT_HOM\tCONTROL_TOTAL_AC\n")
snpfile=open(options.snpfilename, "r")
for line_s1 in snpfile:
	line_s=line_s1.rstrip('\n').split('\t')
	if line_s[0][0]!="#":
		genesnplist=list(set(line_s[1].split(',')))
		sumcounts=sumcount(genesnplist, count_table)
		outfile.write(line_s[0]+"\t"+str(sumcounts[0]-2*sumcounts[1])+"\t"+str(sumcounts[1])+"\t"+str(sumcounts[0])+'\n')
outfile.close()
snpfile.close()

##python count_control.py -v gnomad.test.vcf.gz -o counts.gnomad.txt -s test.out.txt --snpformat CHRPOSREFALT --pop ALL
