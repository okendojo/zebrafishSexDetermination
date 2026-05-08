#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=awking


cd /data/okendojo/zebrafish/data/g3/rna_sequences/meryl_analysis 

awk -v OFS="\t" '{ print $1,$4,$6,$8,$10}' RNA-Seq_only.hapmer_count.txt > kmer_in_reads_and-DBs.txt


