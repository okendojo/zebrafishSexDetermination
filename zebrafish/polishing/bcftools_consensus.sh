#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=bcftoolsConsensus

#load the modules

module load winnowmap
module load samtools
module load meryl
module add minimap2
module add merfin
module add bcftools

WKDIR="/data/okendojo/zebrafish/data/fish11/polish/hybridVariantCall"

cd $WKDIR

bcftools view -Oz fish11_polished_asm.filter.vcf > fish11_polished_asm.filter.vcf.gz

bcftools index fish11_polished_asm.filter.vcf.gz

bcftools consensus -H1 --chain chain_file.txt -f fish11.fasta fish11_polished_asm.filter.vcf.gz > fish11_merfin_polished.fasta
