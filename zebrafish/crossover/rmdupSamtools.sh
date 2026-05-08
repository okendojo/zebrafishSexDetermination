#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=rmdup

cd /data/okendojo/zebrafish/data/G3/TU_AB_alignment

#Load the modules

module load samtools

#See the data paths
samtools rmdup -sS  --output-fmt BAM --verbosity 4 TU_sorted.bam TU_dedup.bam

samtools rmdup -sS  --output-fmt BAM --verbosity 4 AB_sorted.bam AB_dedup.bam








