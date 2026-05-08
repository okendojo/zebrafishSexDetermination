#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=merfinPolish

#load the modules

module load winnowmap
module load samtools
module load meryl
module add minimap2
module add merfin
module add bcftools

cd /data/okendojo/zebrafish/data/fish11/polish/hybridVariantCall

merfin -polish                               \
           -sequence fish11.fasta       \
           -readmers /data/okendojo/zebrafish/data/fish11/polish/asm_polish/hifi_illum_hybrid.meryl   \
	   -loose	\
	   -threads 24	\
           -peak 83.3                        \
           -prob lookup_table.txt            \
           -vcf dv_HYBRID_PACBIO_ILLUMINA/dv_HYBRID_PACBIO_ILLUMINA.vcf.gz  \
           -output fish11_polished_asm
