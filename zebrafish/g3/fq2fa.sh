#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=180:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=fq2fa


#load module

module load seqkit

cd /data/okendojo/zebrafish/data/g3/rna_sequences/untrimmed 

echo "concatenating ont fastq files..............."

cat *.fastq.gz > f2_rna.fastq.gz

echo "Now converting fastq to fasta..............."

#seqkit fq2fa -o f2_rna.fasta -j 26 f2_rna.fastq.gz



