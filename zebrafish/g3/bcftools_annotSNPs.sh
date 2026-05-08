#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=varannot


cd /data/okendojo/zebrafish/data/g3/rna_sequences/haplotype_gene_quantification/g3_variants/haptags_rna_variants/variant_calling

tabix -p vcf 00-All.vcf.gz

# Paths to your input files
input_vcf="merged_top_contigs.vcf.gz"          # Your input VCF file to annotate
dbsnp_vcf="00-All.vcf.gz"       # The dbSNP reference VCF file with rsIDs
output_vcf="annotated_rsID.vcf" # The output VCF with rsIDs

# Annotate the VCF file with rsIDs using bcftools
bcftools annotate -a $dbsnp_vcf -c ID -o $output_vcf -O v $input_vcf


