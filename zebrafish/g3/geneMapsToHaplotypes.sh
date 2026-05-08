#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=200g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=geneMaphaps

#load modules
module add GATK
module load bedtools 


cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants

REF="/data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.dna.primary_assembly.fa"
BED="/data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.111.gtf"

# 1. Filter variants for SNPs

gatk SelectVariants -R ${REF} -V f0_annotatedVariants.vcf.gz -select-type SNP -O snps.vcf

# 2. Filter SNPs for heterozygous calls

gatk SelectVariants -R ${REF} -V snps.vcf -select-type-to-include SNP --restrict-alleles-to BIALLELIC -O het_snps.vcf


# 3. Identify genes mapping to the four haplotypes

bedtools intersect -a het_snps.vcf -b ${BED} -wa -wb > gene_haplotypes_mapping.bed

~                       
