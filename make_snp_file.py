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
parser.add_option("-v", "--vcffile", action="store",dest="vcffilename")
parser.add_option("-o", "--outfile", action="store",dest="outfilename", default="snpfile.txt")
parser.add_option("--genecolname", action="store", dest="genecolname")

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

parser.add_option("--snpformat", action="store",dest="snpformat", default="VCFID")
parser.add_option("--genenull", action="store", dest="genenull", default=".,NA")

options, args = parser.parse_args()

#Try to catch potential errors
if not options.vcffilename:   # if filename is not given
	parser.error('A vcf file is needed')
	sys.exit()

if options.vcffilename.endswith(".gz") is False:   # if vcf filename is not given
	parser.error('Is your vcf file gzipped?')
	sys.exit()

if not options.genecolname:
	parser.error('An INFO field with the gene names to use must be provided')
	sys.exit()

if (options.includevep is not None) or (options.excludevep is not None):
	if not options.vep:
		parser.error('--vep option must be supplied if using VEP annotations')
		sys.exit()
			     
if  options.snpformat!="VCFID" and options.snpformat!="CHRPOSREFALT":   # if filename is not given
	parser.error('SNP format must be "VCFID" or "CHRPOSREFALT"')
	sys.exit()

if options.snponly and options.indelonly:
	parser.error('Please select only --snponly or --indelonly')
	sys.exit()

#Check to make sure all the filters seem well formed
def checkfilter(infofilter):
	if ("[" not in infofilter) or (infofilter.startswith("]")) or (infofilter.endswith("]")) or str(infofilter.split("[")[1].split("]")[0]) not in ["<", ">", "<=", ">=", "=", "!=", "in"]:
		return 0
	else:
		return 1
	
if options.includeinfo is not None:
	for i in range(0, len(options.includeinfo), 1):
		if checkfilter(options.includeinfo[i])==0:
			sys.stdout.write(str(options.includeinfo[i])+" is malformed\n")
			sys.exit()
if options.excludeinfo is not None:
	for i in range(0, len(options.excludeinfo), 1):
		if checkfilter(options.excludeinfo[i])==0:
			sys.stdout.write(str(options.excludeinfo[i])+" is malformed\n")
			sys.exit()
if options.includevep is not None:
	for i in range(0, len(options.includevep), 1):
		if checkfilter(options.includevep[i])==0:
			sys.stdout.write(str(options.includevep[i])+" is malformed\n")
			sys.exit()
if options.excludevep is not None:
	for i in range(0, len(options.excludevep), 1):
		if checkfilter(options.excludevep[i])==0:
			sys.stdout.write(str(options.excludevep[i])+" is malformed\n")
			sys.exit()
		
		 
#Test if something is a number
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def test_include_info(filter, vcfline):
        option_field=filter.split("[")[0]
        option_value=float(filter.split("]")[1])
	if option_field in vcfline:
       		field_value=(";"+vcfline).split((";"+option_field+"="))[1].split(";")[0].split(",")[0]
        	if get_operator_fn(filter.split("[")[1].split("]")[0])(float(field_value), float(option_value)):
			return 1
        	else:
                	return 0
	else:
		return 1

def test_exclude_info(filter, vcfline):
        option_field=filter.split("[")[0]
        option_value=float(filter.split("]")[1])
	if option_field in vcfline:
        	field_value=(";"+vcfline).split((";"+option_field+"="))[1].split(";")[0].split(",")[0]
        	if get_operator_fn(filter.split("[")[1].split("]")[0])(float(field_value), float(option_value)):
			return 0
        	else:
               		return 1
	else:
		return 1

def test_include_vep(filter, vcfline, csq_anno):
	option_field=filter.split("[")[0]
	csq_index=csq_anno.index(option_field)
        option_value=filter.split("]")[1]
	field_value=(";"+vcfline).split((";CSQ="))[1].split(";")[0].split(",")[0].split("|")[csq_index]
	if get_operator_fn(filter.split("[")[1].split("]")[0])(field_value, option_value):
                return 1
      	else:
               	return 0

def test_exclude_vep(filter, vcfline, csq_anno):
        option_field=filter.split("[")[0]
	csq_index=csq_anno.index(option_field)
        option_value=filter.split("]")[1]
	field_value=(";"+vcfline).split(";CSQ=")[1].split(";")[0].split(",")[0].split("|")[csq_index]
	if get_operator_fn(filter.split("[")[1].split("]")[0])(field_value, option_value):
        	return 0
       	else:
               	return 1

def find_vep_gene(genecolname, vcfline, csq_anno):
        csq_index=csq_anno.index(genecolname)
        if "CSQ" in vcfline:
                genename=(";"+vcfline).split(";CSQ=")[1].split(";")[0].split(",")[0].split("|")[csq_index]
        else:
                genename=""
        return genename

def find_info_gene(genecolname, vcfline):
	if genecolname in vcfile:		
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
	vcffile=gzip.open(options.vcffilename, "rb")
	csq_found=0
	for line_vcf1 in vcffile:
		if line_vcf1[0]=="#" and ("ID=CSQ" in line_vcf1):
			csq_anno=line_vcf1.rstrip('\n').replace('"', '').strip('>').split("Format: ")[1].split("|")
			csq_found=1
			break
	if csq_found==0:
		sys.stdout.write("VEP CSQ annotations not found in vcf header\n")
		sys.exit()
	vcffile.close()

if options.vep:
	if options.genecolname not in csq_anno:
		sys.stdout.write("Gene column name not found in VEP annotations\n")
		sys.exit()

##Needs heavy editing still
##Check vcf header to make sure INFO fields are present
##vcffile=gzip.open(options.vcffilename, "rb")	
##	if options.includevep is not None:
##		while iter<len(options.includevep):
	

#Open vcf file
vcffile=BedTool(options.vcffilename)

if options.bedfilename is not None:
	bed=BedTool(options.bedfilename)
	vcffile=vcffile.intersect(bed)

##for line_vcf1 in gzip.open(vcffile.fn):
for line_vcf1 in gzip.open(vcffile.fn):
	line_vcf=line_vcf1.rstrip().split('\t')
	keep=1
	if line_vcf[0][0]!="#":
		if keep==1 and options.passfilter:
			if line_vcf[6]!="PASS":
				keep=0
		if keep==1 and options.snponly:
			if len(line_vcf[3])>1 or len(line_vcf[4])>1:
				keep==0
		if keep==1 and options.indelonly:
			if len(line_vcf[3])==1 and len(line_vcf[4])==1:
				keep==0
  
 #Go through INFO field filters
		if keep==1 and options.includeinfo is not None:
			iter=0
			while keep==1 and iter<len(options.includeinfo):
				filter=options.includeinfo[iter]
				keep=test_include_info(filter, line_vcf[7])
				iter=iter+1
        
		if keep==1 and options.excludeinfo is not None:
			iter=0
			while keep==1 and iter<len(options.excludeinfo):
				filter=options.excludeinfo[iter]
				keep=test_exclude_info(filter, line_vcf[7])
				iter=iter+1
  
   #Go through INFO/VEP field filters
		if keep==1 and options.includevep is not None:
			iter=0
			while keep==1 and iter<len(options.includevep):
				filter=options.includevep[iter]
				keep=test_include_vep(filter, line_vcf[7], csq_anno)
				iter=iter+1
        
		if keep==1 and options.excludevep is not None:
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
	if len(x)>0:
        #Read through hash table and print out variants
        	syn_out=','.join(snptable[x][1])
        	outfile.write(str(x)+"\t"+syn_out+"\n")
outfile.close()

#python make_snp_file.py -o test.out.txt -v gnomad.test.vcf.gz --vep --genecolname SYMBOL --snpformat CHRPOSREFALT --pass --includeinfo "AC[<]5"
