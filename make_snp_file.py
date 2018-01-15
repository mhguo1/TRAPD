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
parser.add_option("--vcffile", action="store",dest="vcffilename")
parser.add_option("--outfile", action="store",dest="outfilename", default="snpfile.txt")

#Filters
parser.add_option("--includeinfo", action="append",dest="includeinfo")
parser.add_option("--excludeinfo", action="append",dest="excludeinfo")
parser.add_option("--includevep", action="append",dest="includevep")
parser.add_option("--excludevep", action="append",dest="excludevep")
parser.add_option("--pass", action="store_true", dest="passfilter")
parser.add_option("--vep", action="store_true", dest="vep")
parser.add_option("--snponly", action="store_true", dest="snponly")
parser.add_option("--indelonly", action="store_true", dest="indelonly")
parser.add_option("--bedfile", action="store", dest="bedfilename")

parser.add_option("--snpformat", action="append",dest="snpformat", default="VCFID")
parser.add_option("--genecolname", action="store", dest="genecolname")
parser.add_option("--snpcolname", action="store", dest="snpcolname", default="SNP")
parser.add_option("--genenull", action="store", dest="genenull", default=".,NA")


options, args = parser.parse_args()

if not options.vcffilename:   # if filename is not given
    parser.error('A vcf file is needed')

if (options.includevep is not None) or (options.excludevep is not None):
	if not options.vep:
		parser.error('--vep option must be supplied if using VEP annotations)

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#Need a better way to handle things in case field is not found in that given line
def test_include_info(filter, vcfline):
        option_field=filter.split("[")[0]
        option_value=float(filter.split("]")[1])
        field_value=(";"+vcfline).split((";"+option_field+"="))[1].split(";")[0].split(",")[0]
        if get_operator_fn(filter.split("[")[1].split("]")[0])(float(field_value), float(option_value)):
		return 1
        else:
                return 0

def test_exclude_info(filter, vcfline):
        option_field=filter.split("[")[0]
        option_value=float(filter.split("]")[1])
        field_value=(";"+vcfline).split((";"+option_field+"="))[1].split(";")[0].split(",")[0]
        if get_operator_fn(filter.split("[")[1].split("]")[0])(float(field_value), float(option_value)):
		return 0
        else:
                return 1

def test_include_vep(filter, vcfline, csq_anno):
	csq_index=csq_anno.index(filter)
        option_field=filter.split("[")[0]
        option_value=filter.split("]")[1]
        field_value=(";"+vcfline).split((";"+option_field+"="))[1].split(";")[0].split(",")[0].split("|")[csq_index]
        if get_operator_fn(filter.split("[")[1].split("]")[0])(float(field_value), float(option_value)):
                return 1
        else:
                return 0


def test_exclude_vep(filter, vcfline, csq_anno):
        csq_index=csq_anno.index(filter)
        option_field=filter.split("[")[0]
        option_value=filter.split("]")[1]
        field_value=(";"+vcfline).split((";"+option_field+"="))[1].split(";")[0].split(",")[0].split("|")[csq_index]
        if get_operator_fn(filter.split("[")[1].split("]")[0])(float(field_value), float(option_value)):
                return 1
        else:
                return 0

def find_vep_gene(genecolname, vcfline, csq_anno):
        csq_index=csq_anno.index(genecolname)
        genename=(";"+vcfline).split(";CSQ=")[1].split(";")[0].split("|")[csq_index]
        return genename

def find_info_gene(genecolname, vcfline):
        genename=(";"+vcfline).split(";"+genecolname+"=")[1].split(";")[0]
        return genename
	
#Function to match operator strings
def get_operator_fn(op):
  return {
	'<' : operator.lt,
	'<=' : operator.le,
	'>' : operator.gt,
	'>=' : operator.gt,
  '=' : operator.eq,
  '!=' : operator.ne,
   }[op]

#Create empty snptable
snptable={}

#Read in vcf header to get VEP CSQ fields
if options.vep:
  if str(options.vcffilename)[-3:]==".gz":
	  vcffile=gzip.open(options.vcffilename, "rb")
  else:
	  vcffile=open(options.vcffilename, "r")
	csq_found=0
	for line_vcf1 in vcffile:
		if line_vcf1[0][0]=="#" and ( "ID=CSQ" in line_vcf1):
			csq_anno=line_vcf1.rstrip('\n').rstrip("\">").split("Format: ")[1].split("|")
			csq_found=1
			break
	if csq_found==0:
		sys.stdout.write("VEP CSQ annotations not found in vcf header\n")
	vcffile.close()



#Open vcf file
vcffile=BedTool(options.vcffilename)

if options.bedfilename is not None:
	bed=BedTool(options.bedfilename)
	vcffile=vcffile.intersect(bed)

for line_vcf1 in open(vcffile.fn):
	line_vcf=line_vcf1.rstrip().split('\t')
	keep=1
	
  if options.passfilter:
	  if line_vcf[6]!="PASS":
			keep=0
	if options.snponly:
		if len(line_vcf[3])>1 or len(line_vcf[4])>1:
			keep==0
	if options.indelonly:
		if len(line_vcf[3])==1 and len(line_vcf[4])==1:
			keep==0
  
 #Go through INFO field filters
	if line_vcf1[0]!="#":
		if options.includeinfo is not None:
			iter=0
			while keep==1 and iter<len(options.includeinfo):
				filter=options.includeinfo[iter]
				keep=test_include_info(filter, line_vcf[7])
				iter=iter+1
        
		if options.excludeinfo is not None:
			iter=0
			while keep==1 and iter<len(options.excludeinfo):
				filter=options.excludeinfo[iter]
				keep=test_exclude_info(filter, line_vcf[7])
				iter=iter+1
  
   #Go through INFO/VEP field filters
		if options.includevep is not None:
			iter=0
			while keep==1 and iter<len(options.includevep):
				filter=options.includevep[iter]
				keep=test_include_vep(filter, line_vcf[7], csq_anno)
				iter=iter+1
        
		if options.excludevep is not None:
			iter=0
			while keep==1 and iter<len(options.excludevep):
				filter=options.excludevep[iter]
				keep=test_exclude_vep(filter, line_vcf[7], csq_anno)
				iter=iter+1
        
#If variant meets all filters, then extract gene name
		if keep==1:
			if options.vep:
				gene=find_vep_gene(options.genecolname, line_vcf[7], csq_anno)
			else:
				gene=find_info_gene(options.genecolname, line_vcf[7])
			if gene not in options.genenull.split(","):
				if options.snpformat=="VCFID":
					snpid=str(line_vcf[2])
				else: 
					snpid=str(line_vcf[0].lstrip("chr"))+":"+str(line_vcf[1])+":"+str(line_vcf[3])+":"+str(line_vcf[4])
				if gene not in snptable:
					snptable[gene]=[gene, []]
				else:
					snptable[gene][1].append(snpid)
pybedtools.cleanup() 			

#Write Output
outfile=open(options.outfilename, "w")
outfile.write("#GENE\tSNPS\n")
for x in snptable:
        #Read through hash table and print out variants
        syn_out=','.join(snptable[x][1])
        outfile.write(str(x)+"\t"+syn_out+"\n")
outfile.close()

#python make_snp_test.py --outfile test.out.txt --vcffile test.vcf.gz --vep --genecolname SYMBOL --pass  --bedfile test.bed
#python make_snp_test.py --outfile test.out.txt --vcffile test.vcf.gz --vep --genecolname SYMBOL --pass --filter "AF[<]0.05" --filter "AC[<]4" --bedfile test.bed
