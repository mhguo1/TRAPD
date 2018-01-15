# PBTest


Finding novel Mendelian disease genes can often be challenging. While some disorders are amenable to family-based analyses (e.g., linkage or segregation), others may be very challenging for a number of reasons. One approach is to use burden testing, where the aggregate burden of rare protein-altering variants in each gene is tested against a set of controls. While one may use a set of available control sequencing data, this is generally too expensive and unavailable in most circumstances. Here, we provide a simple-to-use program that allows for burden testing against publicly available summary level data (e.g., ExAC or gnomAD).

Requirements:
PBTest is written in Python and R. For Python, it is recommended to use Python version 2.7. For R, any version 2.+ should be okay.

Required Python packages:
- optparse
- operator
- pybedtools

Also, if supplying a bed file in Step 1, BEDTools (http://bedtools.readthedocs.io/en/latest/) must also be loaded in the background. 


0) Pre-processing:
There are several pre-processing steps that are necessary before running PBTest: 1) Separating multi-allelic variants, 2) left-aligning indels, 3) annotating your vcf. Below, we provide several sample command lines for performing these steps:


1) Creating a SNP file 
A SNP file is a file that maps qualifying variants to each gene. Qualifying variants are variants that you think may be pathogenic (usually, rare protein-altering variants). You can create a separate SNP file for your cases and controls, or you can make the same SNP file. The SNP file has two columns: 1) Column 1 is your gene name, and 2) Column 2 is a comma separated list of variants assigned to that gene. A variant can be assigned to multiple genes.

To create the SNP file, you must have an annotated vcf from which you define your gene names. The vcf does not need to have any genotypes. 

python make_snp_file.py --vcffile $vcffilename --genecolname $genecol --outfile $outfilename

Required Options
1) --vcffile: This is a path to your vcffile: e.g., /Users/smith/dat/test.vcf.gz. Your vcf does not need to be gzipped, but can be gzipped or bgzipped.

2) --outfile: This is a path to your desired outfile name: e.g., /Users/smith/dat/test.snp_file.txt. 

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
  
 Note that these criteria MUST be surrounded by double quotation marks!

Some examples:
--includeinfo "AC[<]5" #Filters in variants with Allele count (AC) less than five
--excludeinfo "AF[>]0.05" #Filters out variants with allele frequency (AF) greater than 0.05
--includeinfo "SNPEFF_EFFECT[=]missense" #Filters in variants with SNPEFF_EFFECT annotated as missense.

Variants that are kept will meet ALL criteria supplied!

2) VEP INFO field filter options (--includevep, --excludevep). These are criteria based on the VEP INFO field of your vcf to include or exclude variants. They are structured the same as --includeinfo and --excludeinfo, and as many as you want may be used. If these options are used, --vep must be supplied. Some examples:
--includeinfo "BIOTYPE[=]protein_coding" #Include variants where the VEP CSQ consequence is protein_coding
--excludeinfo "consequence[=]synonymous" #Exclude variants where the VEP CSQ consequence is synonymous

3) --snpformat: Format for SNPs. Default is "VCFID". Your SNPs may be defined in any one of two ways.  If you supply the option "VCFID", then the program will use the VCF variant name in column 3 of your vcf (often rsIDs). Alternatively, you may supply "CHRPOSREFALT", in which case variants will be formatted as chr:pos:ref:alt (e.g., 1:1000:A:T).

5) --bedfile: Path to a bed file for regions of interest. Only regions that are inside the bedfile defined regions will be kept. If this option is not supplied, then the entire VCF will be used. Caution that if your chromosome names start in "chr" (e.g., "chr1"), then your bed file should be formatted similarly.


4) --pass: Keep only PASS variants based on the "FILTER" field of your vcf

5) --genenull: Values for which a gene is to be considered null. Default is "NA" or ".". 

6) --vep: Option that should be supplied if you used VEP to annotate your vcf.

7) --snponly: If single nucleotide changes only should be considered.

8) --indelonly: If only indels should be considered.

parser.add_option("--snpcolname", action="store", dest="snpcolname", default="SNP")

parser.add_option("--vep", action="store_true", dest="vep")
parser.add_option("--snponly", action="store_true", dest="snponly")
parser.add_option("--indelonly", action="store_true", dest="indelonly")
parser.add_option("--bedfile", action="store", dest="bedfilename")


