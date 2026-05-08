#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=120g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=bcftools_variantCall

#load module
module load bcftools

cd /data/okendojo/zebrafish/data/g3/rna_sequences/haplotype_gene_quantification/g3_variants

#SET PATHS
Ref_FASTA="/data/okendojo/zebrafish/data/g3/assembly/verkko_asm2/assembly.fasta"


# 1, call variants
bcftools mpileup -Ou --redo-BAQ --min-BQ 20   --threads 24  -f ${Ref_FASTA} AB_T0.markdup.sorted.bam  | bcftools call -mv -Oz --threads 24 -o AB_T0.vcf.gz
