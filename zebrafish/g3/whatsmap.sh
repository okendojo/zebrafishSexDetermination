#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=whatsmap

#load modules
module load samtools
module load whatshap

cd  /data/okendojo/zebrafish/data/g3/rna_sequences/phaser

ref="hapTagsrenamed.fasta"
vcf="haptags_rna_variants/variant_calling/ABTU/ABTU.haplotypecaller.vcf.gz"
BAM="/data/okendojo/zebrafish/data/g3/rna_sequences/phaser/haptags_rna_variants/preprocessing/ABTU/ABTU.markdup.sorted.bam"


#whatshap phase --output ABTU_phased.vcf.gz --sample AB_AB_strain --sample TU_TU_strain --sample TL_TL_strain  --sample WIK_WIK_strain --reference ${ref} ${vcf} $BAM 

whatshap phase --output ABTU_phased.vcf.gz  --recombination-list ab_putativie_recombination.txt --reference ${ref} ${vcf} $BAM 
