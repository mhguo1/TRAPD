# TRAPD


Finding novel Mendelian disease genes can often be challenging. While some disorders are amenable gene discovery to family-based analyses (e.g., linkage or segregation), others may be very challenging for a number of reasons. One approach go identify genes associated with disease is to use burden testing, where the aggregate burden of rare protein-altering variants in each gene is tested against a set of controls. While one may use a set of available control sequencing data, this is generally too expensive and unavailable in most circumstances. Here, we provide a simple-to-use program called TRAPD (Testing Rare vAriants using Public Data) that allows for burden testing against publicly-available summary level data (e.g., ExAC or gnomAD).

Requirements:
TRAPD is written in Python and R. For Python, it is recommended to use Python version 2.7. For R, any version 2.+ should be okay.

Required Python packages:
- optparse
- operator
- pybedtools

BEDTools (http://bedtools.readthedocs.io/en/latest/) must also be loaded in the environment. 


**0) Pre-processing:**
There are several pre-processing steps that are necessary before running TRAPD: 1) Separating multi-allelic variants, 2) left-aligning indels, 3) annotating your vcf. Below, we provide several sample command lines for performing these steps:

	0.1-0.2) Separating multi-allelics and left-aligning indels:
	There are several ways to do this. Please see https://genome.sph.umich.edu/wiki/Variant_Normalization for additional details on this problem. We use Bcftools to accomplish these two steps:
	(https://samtools.github.io/bcftools/bcftools.html):
	bcftools norm -m -any in.vcf.gz | bcftools norm -f Homo_sapiens_assembly19.fasta | bgzip > out.vcf.gz
	
	0.3) Annotation:
	For variant annotation, we have achieved our best results using VEP (https://www.ensembl.org/info/docs/tools/vep/script/index.html). Several additional annotators include SnpEff (http://snpeff.sourceforge.net/) and ANNOVAR (http://annovar.openbioinformatics.org/en/latest/).
	We highly recommend annotating the case and control data in the same way.


**1a) Creating a SNP file**
A SNP file maps qualifying variants to each gene. Qualifying variants are variants that you think may be pathogenic, usually based on minor allele frequency and annotations of predicted effect o the variant (e.g., rare protein-altering variants). You can create a separate SNP file for your cases and controls, or you can make the same SNP file. The SNP file for TRAPD has two columns: 1) Column 1 is your gene name, and 2) Column 2 is a comma separated list of variants assigned to that gene. A variant can be assigned to multiple genes. The header for this file is: "#GENE SNPS".

To create the SNP file, you must have an annotated vcf from which you define your gene names. The vcf does not need to have any genotypes (i.e., can be a sites-only vcf). 

python make_snp_file.py --vcffile $vcffilename --genecolname $genecol --outfile $outfilename

Required Options
1) -v, --vcffile: This is a path to your vcffile: e.g., /Users/smith/dat/test.vcf.gz. Your vcf must be gzipped or bgzipped.

2) -o, --outfile: This is a path to your desired outfile name: e.g., /Users/smith/dat/snp_file.txt. The default is "snpfile.txt"

3) --genecolname: This is a field within the INFO field of your vcf telling the script which gene that variant belongs to. For SNPEFF, this is typically SNPEFF_GENENAME. If you used VEP to annotate your vcf (see Step 0), you must supply the --vep option below, and you'll use the column name within the CSQ field for VEP (usually "genename" or the like).

Additional Options:
1) INFO field filter options (--includeinfo, --excludeinfo): These are criteria based on the INFO field of your vcf to include or exclude variants. You may include as many inclusion or exclusion criteria as you wish. 

These criteria are structured as: "FIELD[operator]threshold". Here, FIELD is any field within the INFO field of your vcf (e.g., AC, ExAC_AF, AF). Operator is any operator from the following: 
	'<': less than
  	'<=' : less than or equal to
	'>' : greater than
	'>=' : greater than or equal to
 	 '=' : equals
  	'!=' : does not equal
	'in': in
  
 Also, note that if you use the "in" operator, you should include a list of test values enclosed by parentheses (e.g., "annotation[in](missense,nonsense,frameshift)"

Some examples:
--includeinfo "AC[<]5" #Filters in variants with Aalele count (AC) less than five
--excludeinfo "AF[>]0.05" #Filters out variants with allele frequency (AF) greater than 0.05
--includeinfo "SNPEFF_EFFECT[=]missense" #Filters in variants with SNPEFF_EFFECT annotated as missense.

Variants that are kept will meet ALL criteria supplied!

2) VEP INFO field filter options (--includevep, --excludevep). These are criteria based on the VEP INFO field of your vcf to include or exclude variants. They are structured the same as --includeinfo and --excludeinfo, and as many as you want may be used. If these options are used, --vep must be supplied. Some examples:
--includeinfo "BIOTYPE[=]protein_coding" #Include variants where the VEP CSQ consequence is protein_coding
--excludeinfo "consequence[=]synonymous" #Exclude variants where the VEP CSQ consequence is synonymous

3) --snpformat: Format for SNPs. Default is "VCFID". Your SNPs may be defined in any one of two ways.  If you supply the option "VCFID", then the program will use the VCF variant name in column 3 of your vcf (often rsIDs). Alternatively, you may supply "CHRPOSREFALT", in which case variants will be formatted as chr:pos:ref:alt (e.g., 1:1000:A:T).

4) --bedfile: Path to a bed file for regions of interest. Only regions that are within the bed file-defined regions will be kept. If this option is not supplied, then the entire VCF will be used. Caution that if your chromosome names start in "chr" (e.g., "chr1"), then your bed file should be formatted similarly.

5) --pass: Keep only PASS variants based on the "FILTER" field of your vcf

6) --genenull: Values for which a gene is to be considered null. Default is "NA" or ".". Genes names having any of these values will be excluded from teh output.

7) --vep: Option that should be supplied if you used VEP to annotate your vcf.

8) --snponly: If single nucleotide changes only should be considered.

9) --indelonly: If only indels should be considered.

Output: The output file will contain two columns: 1) Column 1 will be a list of genes and will have the header "#GENE", and 2) Column 2 will be a comma separated list of SNPs assigned to that gene with the header "SNPS".

**1b) Merging SNP files**
This is an optional step if you need to merge two SNP files (for example, if you performed step 1a separately for SNPs and indels). It can also be used if you perfomed Step 1a separately for each chromosome. It has two required options:

Required Options
1) -s, --snpfiles: This is a comma-separated list of SNP files from Step 1a you are trying to merge.

2) -o, --outfile: This is a path to your desired outfile name: e.g., /Users/smith/dat/out.txt.


**2a) Counting carriers in case cohort**
This script will tabulate the number of cases carrying qualifying variants in each gene as defined by a SNP file. The script will generate three counts for each gene:
- CASE_COUNT_HET: # of individuals carrying one (and exactly one) heterozygous qualifying variant in the gene
- CASE_COUNT_CH: # of individuals carrying at least two heterozygous qualifying variants in the gene
- CASE_COUNT_HOM: # of individuals carrying at least one homozygous qualifying variant in the gene.

The command takes in a vcf file containing case sample genotypes and a SNP file listing the qualifying variants for each gene. The general command is:
python count_cases.py -v test.vcf.gz -s snpfile.txt -o controlcounts.txt [--snpformat --samplefile --pass --maxAC --maxAF --GTfield]. 

Required Options
1) -v, --vcffile: This is the path to your VCF file containing case genotypes: e.g., /Users/smith/dat/test.vcf.gz. Your vcf must be gzipped or bgzipped.

2) -s, --snpfile: This is the path to your SNP file containing mappings of qualifying variants to gene (see Step 1). 

3) -o, --outfile: This is a path to your desired outfile name: e.g., /Users/smith/dat/out.txt. The default is "case_counts.txt"

Additional Options:
4) --snpformat: Format for SNPs. Default is "VCFID". Your SNPs may be defined in any one of two ways.  If you supply the option "VCFID", then the program will use the VCF variant name in column 3 of your vcf (often rsIDs). Alternatively, you may supply "CHRPOSREFALT", in which case variants will be formatted as chr:pos:ref:alt (e.g., 1:1000:A:T).

5) --samplefile: Optional file containing list of samples to use. The file should contain one sample per row. Only samples in this list and in the VCF will be used. 

6) --pass: Keep only PASS variants based on the "FILTER" field of your vcf

7) --maxAC: Keep only variants with allele count (AC) less than this value. Note that this is calculated based on all samples in the VCF. The default is 99999.

8) --maxAF: Keep only variants with allele frequency (AF) less than this value. Note that this is calculated based on all samples in the VCF. The default is 1.0.

9) --minAN: Keep only variants with allele number (AN) greater than this value. The default is 0. 

10) --popmaxAF: Keep only variants with population maximum allele frequency (AF) less than this value. If gnomAD or ExAC is used, this option will use the values from the downloaded vcf. If "generic" database is used, there must be a field in "INFO" called "AF_POPMAX". The default is 1.0.

11) --GTfield: The format field within the genotype data from which genotypes should be extracted. The default is "GT"

11) --bedfile: Path to a bed file for regions of interest. Only regions that are within the bed file-defined regions will be kept. If this option is not supplied, then the entire VCF will be used. Caution that if your chromosome names start in "chr" (e.g., "chr1"), then your bed file should be formatted similarly.

Output: The output file will contain three columns: Column 1 will be a list of genes and will have the header "#GENE", Column 2 will be the number of heterozygous cases with the header "CASE_COUNT_HET", Column 3 will be the number of potentially compound heterozygous cases (those carrying 2 variants in the gene), with the header "CASE_COUNT_CH", and Column 4 will be the number of homozygous cases with the header "CASE_COUNT_HOM", and column 5 will be the total AC in the gene with the header "CASE_TOTAL_AC".

**2b) Counting carriers in public control cohorts**
This script will tabulate the approximate number of controls carrying qualifying variants in each gene as defined by a SNP file. Currently, this script has been configured to run using ExAC (http://exac.broadinstitute.org/downloads) or gnomAD (http://gnomad.broadinstitute.org/) data. The script will generate two counts for each gene:
- CONTROL_COUNT_HET: Sum of allele counts of heterozygous qualifying variants in a given gene. 
- CONTROL_COUNT_HOM: Sum of individuals carrying homozygous qualifying variants in a given gene. 


Required Options
1) -v, --vcffile: This is the path to  VCF file containing control genotypes: e.g., /Users/smith/dat/public.vcf.gz. The vcf must be gzipped or bgzipped. The VCF should be downloaded from ExAC (http://exac.broadinstitute.org/downloads) or gnomAD (http://gnomad.broadinstitute.org/)

2) -s, --snpfile: This is the path to your SNP file containing mappings of qualifying variants to gene (see Step 1). 

3) -o, --outfile: This is a path to your desired outfile name: e.g., /Users/smith/dat/out.txt. The default is "case_counts.txt"

Additional Options:
4) --snpformat: Format for SNPs. Default is "VCFID". Your SNPs may be defined in any one of two ways.  If you supply the option "VCFID", then the program will use the VCF variant name in column 3 of your vcf (often rsIDs). Alternatively, you may supply "CHRPOSREFALT", in which case variants will be formatted as chr:pos:ref:alt (e.g., 1:1000:A:T).

5) --pop: Comma separated list of continental populations to use. For ExAC, these include AFR, AMR, EAS, FIN, NFE, SAS, OTH.  For gnomad, these include AFR, AMR, ASJ, EAS, FIN, NFE, SAS, OTH. If ALL is included, then all populations are used. The default is "ALL"

6) -d, --database: Control database used. Default is "generic" which is any vcf that contains at least AC and AN in the INFO column. Other options are "gnomad" or "exac"

7) --pass: Keep only PASS variants based on the "FILTER" field of VCF

8) --maxAC: Keep only variants with allele count (AC) less than this value. Note that this is based on the INFO/AC field in the VCF. The default is 99999.

9) --maxAF: Keep only variants with allele frequency (AF) less than this value. Note that this is based on the INFO/AF field in the VCF. The default is 1.0.

10) --minAN: Keep only variants with allele number (AN) greater than this value.  Note that this is based on the INFO/AF field in the VCF. The default is 0. 

11) --homcol: This argument is only relevant if the database is "generic". This argument specifies the INFO field that contains the counts of the number of homozygotes. If not specified for a generic database, then homozygotes will not be counted.

12) --bedfile: Path to a bed file for regions of interest. Only regions that are within the bed file-defined regions will be kept. If this option is not supplied, then the entire VCF will be used. Caution that if your chromosome names start in "chr" (e.g., "chr1"), then your bed file should be formatted similarly.

Output: The output file will contain three columns: Column 1 will be a list of genes and will have the header "#GENE", Column 2 will be the number of heterozygous controls, Column 3 will be the number of homozygous controls, and Column 4 will be the total AC 


**3) Run burden testing**
This script will run the actual burden testing. It performs a one-sided Fisher's exact test to determine if there is a greater burden of qualifying variants in cases as compared to controls for each gene. It will perform this burden testing under a dominant and a recessive model.

It requires R; the script was tested using R v3.1, but any version of R should work.

The script has 5 required options:
1) --casefile: Path to the counts file for the cases, as generated in Step 2A
2) --casesize: Number of cases that were tested in Step 2A
3) --controlfile: Path to the counts file for the controls, as generated in Step 2B
4) --controlsize: Number of controls that were tested in Step 2B. If using ExAC or gnomAD, please refer to the respective documentation for total sample size
5) --output: Output file path/name

Output: A tab delimited file with 10 columns: Column 1 ("#GENE") will be a list of genes; Column 2 ("CASE_COUNT_HET") will be the number of heterozygous cases; Column 3 ("CASE_COUNT_CH") will be the number of potentially compound heterozygous cases (those carrying 2 variants in the gene); Column 4 ("CASE_COUNT_HOM") will be the number of homozygous cases; Column 5 ("CASE_TOTAL_AC") will be the total AC in the gene for cases; Column 6 will be the estimated number of heterozygous controls ("CONTROL_COUNT_HET"); Column 7 will be the number of homozygous controls ("CONTROL_COUNT_HOM"), and Column 8 will be the total AC in controls ("CONTROL_TOTAL_AC"), Column 9 ("P_DOM") will be the p-value under the dominant model, and Column 10 ("P_REC") will be the p-value under the recessive model.



Citing TRAPD:
Guo MH, Plummer L, Chan Y-M, Hirschhorn JN, Lippincott MF. Burden testing of rare variants identified through exome sequencing using publicly available control data. American Journal of Human Genetics. 2018 (manuscript in press).


