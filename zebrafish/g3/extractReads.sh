#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=1000g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=extract

cd /data/okendojo/zebrafish/data/g3/rna_sequences/haplotype_gene_quantification

#grep "AB" forwardTags.csv > ab_forward_reads.csv
#grep "TU" forwardTags.csv > tu_forward_reads.csv
#grep "TL" forwardTags.csv > tl_forward_reads.csv
#grep "WIK" forwardTags.csv > wik_forward_reads.csv

#grep "AB" reverseTags.csv > ab_reverse_reads.csv
#grep "TU" reverseTags.csv  > tu_reverse_reads.csv 
#grep "TL" reverseTags.csv  > tl_reverse_reads.csv
#grep "WIK" reverseTags.csv > wik_reverse_reads.csv


#cut -f1 tmp.txt | sort | uniq -c > strain_readsnumber.csv

for fastq in *.fastq ; do gzip $fastq ; done 
