#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=90:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=indexing

#load modules
module load bowtie

cd /data/okendojo/zebrafish/data/g3/mir_430Gene/cage_data/asm_indices

bowtie2 -p 24 -q --local -x /data/okendojo/zebrafish/data/g3/mir_430Gene/cage_data/asm_indices/fish11_index -U /data/okendojo/zebrafish/data/g3/mir_430Gene/cage_data/fastq/SRX156356_SRR516560.fastq.gz -b /data/okendojo/zebrafish/data/g3/mir_430Gene/cage_data/fish11/SRX156356_SRR516560_fish11.bam 

bowtie2 -p 24 -q --local -x /data/okendojo/zebrafish/data/g3/mir_430Gene/cage_data/asm_indices/fish11_index -U /data/okendojo/zebrafish/data/g3/mir_430Gene/cage_data/fastq/SRX156356_SRR516561.fastq.gz -b /data/okendojo/zebrafish/data/g3/mir_430Gene/cage_data/fish11/SRX156356_SRR516561_fish11.bam

