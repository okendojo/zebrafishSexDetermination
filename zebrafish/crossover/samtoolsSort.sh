#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=samtoolsRun

cd /data/okendojo/zebrafish/data/G3/TU_AB_alignment

#Load the modules

module load samtools

#See the data paths
REF_GENOME=/data/okendojo/zebrafish/refGenome/GRCz11_genomic.fasta

#Convert sam to bam
#samtools view -b -u -o TU.bam -t ${REF_GENOME}.fai TU.sam
#samtools view -b -u -o AB.bam -t ${REF_GENOME}.fai AB.sam

#Sort bam files
samtools sort -u -m 2G -o AB_sorted.bam  --threads 20 -O BAM AB.bam
samtools sort -u -m 2G -o TU_sorted.bam  --threads 20 -O BAM TU.bam











