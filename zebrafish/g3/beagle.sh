#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=beaglePhasing

#Load modules

module load Beagle
module load bcftools

cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants

bcftools view -s "AB_AB_strain" f0_annotatedVariants.vcf.gz -o AB.vcf.gz -Oz 

bcftools view -s "TU_TU_strain" f0_annotatedVariants.vcf.gz -o TU.vcf.gz -Oz

bcftools view -s "TL_TL_strain" f0_annotatedVariants.vcf.gz -o TL.vcf.gz -Oz

bcftools view -s "WIK_WIK_strain" f0_annotatedVariants.vcf.gz -o WIK.vcf.gz -Oz

#Identify how many SNPs each haplotype has 

# 1. Phase the genome using beagle
##java -jar $BEAGLE_JAR gt=f0_annotatedVariants.vcf.gz out=phasedParentals.out nthreads=24

# 2. Extract Haplotype Information:
##bcftools query -i 'ALT="SNP" && TYPE="snp"' -f '%CHROM\t%POS\t%HAP\t%REF\t%ALT\n' your_phased.vcf > haplotype_info.txt

# 3. Count SNPs for Each Haplotype and Gene:
##import pandas as pd

# Load haplotype information into a DataFrame
##hap_df = pd.read_csv('haplotype_info.txt', sep='\t', header=None, names=['CHROM', 'POS', 'HAP', 'REF', 'ALT'])

# Assuming you have a separate file with gene information (gene_bed_file.bed)
##gene_bed = pd.read_csv('gene_bed_file.bed', sep='\t', header=None, names=['CHROM', 'START', 'END', 'GENE'])

# Merge haplotype and gene information based on chromosome
###merged_df = pd.merge(hap_df, gene_bed, on='CHROM')

# Filter SNPs within gene boundaries
##gene_snps = merged_df[(merged_df['POS'] >= merged_df['START']) & (merged_df['POS'] <= merged_df['END'])]

# Count SNPs per haplotype and gene
##snp_counts = gene_snps.groupby(['GENE', 'HAP']).size().reset_index(name='SNP_COUNT')

##print(snp_counts)

