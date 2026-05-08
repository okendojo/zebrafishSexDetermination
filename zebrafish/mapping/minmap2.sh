#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=map

#Load the module
module load minimap2 samtools

cd /data/okendojo/zebrafish/data/assembliesVariants

minimap2 -t 30 -ax map-ont /data/okendojo/zebrafish/data/fish11/finalPolish/nhgri_polished_fish11.fasta  fish11_rDNA.fasta | samtools view -Sb - | samtools sort -o F11_rDNA_sorted_aligned.bam

minimap2 -t 30 -ax map-ont /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/asm_polished/NHGRI_Fish6_cons.fasta  f6_rna.fa | samtools view -Sb - | samtools sort -o F6_rDNA_sorted_aligned.bam
