#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=bamTobed



#Load the module
module load deeptools
module load samtools
module load bedtools


cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/ncb_ref_analysis


for file in *.bam ; do bedtools bamtobed -i ${file} > ${file}.bed ; done 


