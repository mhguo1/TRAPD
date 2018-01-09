#!/usr/bin/python
import optparse
import sys
import gzip 

#Parse options
parser = optparse.OptionParser()
parser.add_option("-s", "--snpfile", action="store",dest="snpfilename")
parser.add_option("-v", "--vcffile", action="store",dest="vcffilename")
parser.add_option("-o", "--outfile", action="store",dest="outfilename", default="out.txt") 

parser.add_option("--snpformat", action="store",dest="snpformat", default="VCFID")
parser.add_option("--snpcolname", action="store",dest="snpcolname", default="NA")
parser.add_option("--samplefile", action="store",dest="samplefilename", default="ALL")
parser.add_option("--recessive", action="store_true",dest="recessive")
parser.add_option("--maxAF", action="store",dest="maxAF", default=1)
parser.add_option("--maxAC", action="store",dest="maxAC", default=99999)
parser.add_option("--GTfield", action="store",dest="gtfield", default="GT")

options, args = parser.parse_args()
if not options.snpfilename:   # if filename is not given
    parser.error('A file containing a list of SNPs is needed')

if not options.vcffilename:   # if vcf filename is not given
    parser.error('A vcf file is needed')

def findhetcarriers(vcfline, gtname, snplist, snpformat):
	#Find the column in the genotype field corresponding to the genotypes
	gtcol=vcfline.split('\t')[8].split(":").index(gtname)

	if snpformat=="VCFID":
		snpid=vcfline.split('\t')[2]
	else:
		snpid=str(vcfline.split('\t')[0]).lstrip("chr")+":"+str(vcfline.split('\t')[1])+":"+str(vcfline.split('\t')[3])+":"+str(vcfline.split('\t')[4])
	
	#Extract genotypes 
	gt=[i.split(':')[gtcol] for i in vcfline.rstrip().split('\t')[9:]]

	#Find carriers
	carriers=[i for i,val in enumerate(gt) if str(val) in ["0/1", "1/0", "1/1","0|1", "1|0", "1|1"]]

	return carriers

def findhomcarriers(vcfline, gtname, snplist, snpformat):
	#Find the column in the genotype field corresponding to the genotypes
	gtcol=vcfline.split('\t')[8].split(":").index(gtname)

	if snpformat=="VCFID":
		snpid=vcfline.split('\t')[2]
	else:
		snpid=str(vcfline.split('\t')[0]).lstrip("chr")+":"+str(vcfline.split('\t')[1])+":"+str(vcfline.split('\t')[3])+":"+str(vcfline.split('\t')[4])
	
	#Extract genotypes 
	gt=[i.split(':')[gtcol] for i in vcfline.rstrip().split('\t')[9:]]

	#Find carriers
	hetcarriers=[i for i,val in enumerate(gt) if str(val) in ["0/1", "1/0", "0|1", "1|0"]]
	homcarriers=[i for i,val in enumerate(gt) if str(val) in ["0/1", "1/0","1/1", "0|1", "1|0", "1|1"]]
	carriers=hetcarriers+homcarriers+homcarriers
	return carriers

def findsampleindex(vcfline, samplefilename):
	#This takes the vcf header line and finds the indices corresponding to the individuals present in the sample file
	samplenames=vcfline.rstrip().split('\t')[9:]

	#If User doesn't provide sample list, assume all samples in vcf
	if samplefilename=="ALL":
		sampleindex=range(0, len(samplenames),1)

	#else, find the indices corresponding to the samples in the user-provided list
	else:
		#Generate sample list
		sample_list=[]
		sample_file=open(samplefilename, "r")
		for line_s1 in sample_file:
        			sample_list.append(line_s1.rstrip())
		sample_file.close()
		sampleindex=[i for i,val in enumerate(samplenames) if str(val) in samplelist]
	return sampleindex

def makesnplist(snpfile, snpcolname):
	#Makes a list of SNPs present in the snpfile
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


def calculatecount(genesnps, snptable):
	#This will generate an aggregate count for a given gene.
        gt_index=[]
        for s in range(0, len(genesnps), 1):
                if genesnps[s] in snptable:
                        tempsnp=genesnps[s]
                        if type(tempsnp) is list:
                                gt_index=gt_index+snptable[tempsnp][1]
                        else:
                                gt_index.append(tempsnp)
        ac=len(list(set(gt_index)))
        return ac


listofsnps=makesnplist(options.snpfilename, options.snpcolname)
tableout={}
#Open vcf file
if str(options.vcffilename)[-3:]==".gz":
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
			if options.recessive:
				templist=findhomcarriers(line_vcf1, options.gtfield, listofsnps, options.snpcolname)
				tableout[snpid]=[snpid, len(list(set(sampleindices) & set([x for x in templist if templist.count(x) > 1])))]
			else:
				templist=findhetcarriers(line_vcf1, options.gtfield, listofsnps, options.snpcolname)
				tableout[snpid]=[snpid, list(set(sampleindices) & set(templist))]
	
	#Find indices of samples in the sample file
	elif line_vcf[0]=="#CHROM":
		sampleindices=findsampleindex(line_vcf1, options.samplefilename)

vcffile.close()


#Test by writing output
outfile=open(options.outfilename, "w")
snpfile=open(options.snpfilename, "r")
for line_s1 in snpfile:
	line_s=line_s1.rstrip('\n').split('\t')
	if line_s[0]!="GENE":
		#snplist=list(set(line_s[1].split(',')))
		snplist=line_s[1].split(',')
		count=calculatecount(snplist, tableout)
		outfile.write(line_s[0]+"\t"+str(count)+'\n')
outfile.close()
snpfile.close()


