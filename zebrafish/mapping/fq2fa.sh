#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=27:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=fq2fa


#Load modulels
module load seqkit

cd /data/okendojo/zebrafish/data/g3/sex_project

seqkit fq2fa -o 10ww.fasta -j 24 10WW/m84270_240823_162003_s2.hifi_reads.fastq.gz

seqkit fq2fa -o 1zz.fasta -j 24 1ZZ/m84270_240821_131130_s1.hifi_reads.fastq.gz 

seqkit fq2fa -o 2zz.fasta -j 24 2ZZ/m84270_240821_154547_s2.hifi_reads.fastq.gz 

seqkit fq2fa -o 4zz.fasta -j 24 4ZZ/m84270_240821_205429_s4.hifi_reads.fastq.gz

seqkit fq2fa -o 6ww.fasta -j 24 6WW/m84270_240823_134527_s1.hifi_reads.fastq.gz
