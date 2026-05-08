#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=LDheat

#Load the modules
module load R

cd /data/okendojo/zebrafish/data/g3/rna_sequences/haplotype_gene_quantification/g3_variants/haptags_rna_variants/variant_calling

Rscript ldBlock.R
