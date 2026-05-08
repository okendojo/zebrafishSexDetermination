#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=contigSize

cd /data/okendojo/zebrafish/data/g3/illumina/meryldbs

cut -f5 g3_asm_Haploblocks_2.bed | sort | uniq -c  > contig_size_counts.txt 
