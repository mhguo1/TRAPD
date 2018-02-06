# TRAPD


Finding novel Mendelian disease genes can often be challenging. While some disorders are amenable to family-based analyses (e.g., linkage or segregation), others may be very challenging for a number of reasons. One approach is to use burden testing, where the aggregate burden of rare protein-altering variants in each gene is tested against a set of controls. While one may use a set of available control sequencing data, this is generally too expensive and unavailable in most circumstances. Here, we provide a simple-to-use program called TRAPD (Testing Rare vAriants using Public Data) that allows for burden testing against publicly available summary level data (e.g., ExAC or gnomAD).

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
	There are several ways to do this. Please see https://genome.sph.umich.edu/wiki/Variant_Normalization for additional details on this problem. We use Bcftools to accomplish these two steps using Bcftools (https://samtools.github.io/bcftools/bcftools.html):
	bcftools norm -m -any in.vcf.gz | bcftools norm -f Homo_sapiens_assembly19.fasta | bgzip > out.vcf.gz
	
	0.3) Annotation:
	For variant annotation, we typically use VEP (https://www.ensembl.org/info/docs/tools/vep/script/index.html) and have achieved the best results with VEP. Several additional annotators include SnpEff (http://snpeff.sourceforge.net/) and ANNOVAR (http://annovar.openbioinformatics.org/en/latest/).
	We highly recommend annotating the case and control data in the same way.


**1) Creating a SNP file**
A SNP file maps qualifying variants to each gene. Qualifying variants are variants that you think may be pathogenic (usually, rare protein-altering variants). You can create a separate SNP file for your cases and controls, or you can make the same SNP file. The SNP file has two columns: 1) Column 1 is your gene name, and 2) Column 2 is a comma separated list of variants assigned to that gene. A variant can be assigned to multiple genes. The header for this file is: "#GENE SNPS".

To create the SNP file, you must have an annotated vcf from which you define your gene names. The vcf does not need to have any genotypes. 

python make_snp_file.py --vcffile $vcffilename --genecolname $genecol --outfile $outfilename

Required Options
1) -v, --vcffile: This is a path to your vcffile: e.g., /Users/smith/dat/test.vcf.gz. Your vcf must be gzipped or bgzipped.

2) -o, --outfile: This is a path to your desired outfile name: e.g., /Users/smith/dat/snp_file.txt. The default is "snpfile.txt"

3) --genecolname: This is a field within the INFO field of your vcf telling the script which gene that variant belongs to. For SNPEFF, this is typically SNPEFF_GENENAME. If you used VEP to annotate your vcf (see Step 0), you must supply the --vep option below, and you'll use the column name within the CSQ field for VEP (usually "genename" or the like).

Additional Options
1) INFO field filter options (--includeinfo, --excludeinfo): These are criteria based on the INFO field of your vcf to include or exclude variants. You may include as many inclusion or exclusion criteria as you wish. 

These criteria are structured as: "FIELD[operator]threshold". Here, FIELD is any field within the INFO field of your vcf (e.g., AC, ExAC_AF, AF). Operator is any operator in: 
	'<': less than
  	'<=' : less than or equal to
	'>' : greater than
	'>=' : greater than or equal to
 	 '=' : equals
  	'!=' : does not equal
	"in": in
  
 Note that these criteria MUST be surrounded by double quotation marks!
 Also, note that if you use the "in" operator, you should include a list of test values enclosed by parentheses (e.g., "annotation[in](missense,nonsense,frameshift)"

Some examples:
--includeinfo "AC[<]5" #Filters in variants with Allele count (AC) less than five
--excludeinfo "AF[>]0.05" #Filters out variants with allele frequency (AF) greater than 0.05
--includeinfo "SNPEFF_EFFECT[=]missense" #Filters in variants with SNPEFF_EFFECT annotated as missense.

Variants that are kept will meet ALL criteria supplied!

2) VEP INFO field filter options (--includevep, --excludevep). These are criteria based on the VEP INFO field of your vcf to include or exclude variants. They are structured the same as --includeinfo and --excludeinfo, and as many as you want may be used. If these options are used, --vep must be supplied. Some examples:
--includeinfo "BIOTYPE[=]protein_coding" #Include variants where the VEP CSQ consequence is protein_coding
--excludeinfo "consequence[=]synonymous" #Exclude variants where the VEP CSQ consequence is synonymous

3) --snpformat: Format for SNPs. Default is "VCFID". Your SNPs may be defined in any one of two ways.  If you supply the option "VCFID", then the program will use the VCF variant name in column 3 of your vcf (often rsIDs). Alternatively, you may supply "CHRPOSREFALT", in which case variants will be formatted as chr:pos:ref:alt (e.g., 1:1000:A:T).

5) --bedfile: Path to a bed file for regions of interest. Only regions that are within the bed file-defined regions will be kept. If this option is not supplied, then the entire VCF will be used. Caution that if your chromosome names start in "chr" (e.g., "chr1"), then your bed file should be formatted similarly.

4) --pass: Keep only PASS variants based on the "FILTER" field of your vcf

5) --genenull: Values for which a gene is to be considered null. Default is "NA" or ".". 

6) --vep: Option that should be supplied if you used VEP to annotate your vcf.

7) --snponly: If single nucleotide changes only should be considered.

8) --indelonly: If only indels should be considered.

Output: The output file will contain two columns: 1) Column 1 will be a list of genes and will have the header "#GENE", and 2) Column 2 will be a comma separated list of SNPs assigned to that gene with the header "SNPS".



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

Additional Options
4) --snpformat: Format for SNPs. Default is "VCFID". Your SNPs may be defined in any one of two ways.  If you supply the option "VCFID", then the program will use the VCF variant name in column 3 of your vcf (often rsIDs). Alternatively, you may supply "CHRPOSREFALT", in which case variants will be formatted as chr:pos:ref:alt (e.g., 1:1000:A:T).

5) --samplefile: Optional file containing list of samples to use. File should contain one sample per row. Only samples in this list and in the VCF will be used. 

6) --pass: Keep only PASS variants based on the "FILTER" field of your vcf

7) --maxAC: Keep only variants with allele count (AC) less than this value. Note that this is calculated based on all samples in the VCF (i.e., the INFO field is not used). The default is 99999.

8) --maxAF: Keep only variants with allele frequency (AF) less than this value. Note that this is calculated based on all samples in the VCF (i.e., the INFO field is not used). The default is 1.0.

9) --GTfield: The format field within the genotype data from which genotypes should be extracted. The default is "GT"




**2b) Counting carriers in public control cohorts**
This script will tabulate the approximate number of controls carrying qualifying variants in each gene as defined by a SNP file. Currently, this script has been configured to run using ExAC (http://exac.broadinstitute.org/downloads) or gnomAD (http://gnomad.broadinstitute.org/) data. The script will generate two counts for each gene:
- CONTROL_COUNT_HET: Sum of allele counts of heterozygous qualifying variants in a given gene. 
- CONTROL_COUNT_HOM: Sum of individuals carrying homozygous qualifying variants in a given gene. 


Required Options
1) -v, --vcffile: This is the path to  VCF file containing control genotypes: e.g., /Users/smith/dat/public.vcf.gz. The vcf must be gzipped or bgzipped. The VCF should be downloaded from ExAC (http://exac.broadinstitute.org/downloads) or gnomAD (http://gnomad.broadinstitute.org/)

2) -s, --snpfile: This is the path to your SNP file containing mappings of qualifying variants to gene (see Step 1). 

3) -o, --outfile: This is a path to your desired outfile name: e.g., /Users/smith/dat/out.txt. The default is "case_counts.txt"

Additional Options
4) --snpformat: Format for SNPs. Default is "VCFID". Your SNPs may be defined in any one of two ways.  If you supply the option "VCFID", then the program will use the VCF variant name in column 3 of your vcf (often rsIDs). Alternatively, you may supply "CHRPOSREFALT", in which case variants will be formatted as chr:pos:ref:alt (e.g., 1:1000:A:T).

5) --pop: Comma separated list of continental populations to use. For ExAC, these include AFR, AMR, EAS, FIN, NFE, SAS, OTH.  For gnomad, these include AFR, AMR, ASJ, EAS, FIN, NFE, SAS, OTH. If ALL is included, then all populations are used. The default is "ALL"

6) --pass: Keep only PASS variants based on the "FILTER" field of VCF

7) --maxAC: Keep only variants with allele count (AC) less than this value. Note that this is based on the INFO/AC field in the VCF. The default is 99999.

8) --maxAF: Keep only variants with allele frequency (AF) less than this value. Note that this is based on the INFO/AF field in the VCF. The default is 1.0.
