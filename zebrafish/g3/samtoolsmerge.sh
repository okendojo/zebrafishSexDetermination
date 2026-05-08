#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=200g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=samtoolsMerge

#load modules
module add bedtools
module load samtools

cd /data/okendojo/zebrafish/data/g3/rna_sequences/haplotype_gene_quantification/g3_variants

samtools merge -@ 24 -O BAM -o g3_merged.bam haptags_rna_variants/preprocessing/10TL_T6/10TL_T6.markdup.sorted.bam haptags_rna_variants/preprocessing/11TL_T12/11TL_T12.markdup.sorted.bam haptags_rna_variants/preprocessing/12TL_T24/12TL_T24.markdup.sorted.bam haptags_rna_variants/preprocessing/13WIK_T0/13WIK_T0.markdup.sorted.bam  haptags_rna_variants/preprocessing/14WIK_T6/14WIK_T6.markdup.sorted.bam haptags_rna_variants/preprocessing/15WIK_T12/15WIK_T12.markdup.sorted.bam haptags_rna_variants/preprocessing/16WIK_T24/16WIK_T24.markdup.sorted.bam haptags_rna_variants/preprocessing/1AB_T0/1AB_T0.markdup.sorted.bam haptags_rna_variants/preprocessing/2AB_T6/2AB_T6.markdup.sorted.bam haptags_rna_variants/preprocessing/3AB_T12/3AB_T12.markdup.sorted.bam haptags_rna_variants/preprocessing/4AB_T24/4AB_T24.markdup.sorted.bam haptags_rna_variants/preprocessing/5TU_T0/5TU_T0.markdup.sorted.bam haptags_rna_variants/preprocessing/6TU_T6/6TU_T6.markdup.sorted.bam haptags_rna_variants/preprocessing/7TU_T12/7TU_T12.markdup.sorted.bam haptags_rna_variants/preprocessing/8TU_T244/8TU_T244.markdup.sorted.bam haptags_rna_variants/preprocessing/9TL_T0/9TL_T0.markdup.sorted.bam 


# Index the merged bam
samtools index g3_merged.bam

