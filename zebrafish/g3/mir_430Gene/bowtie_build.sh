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

bowtie2-build /data/okendojo/zebrafish/data/g3/mir_430Gene/cage_data/assembly/chr4.fasta chr4_index 

#bowtie-build /data/okendojo/zebrafish/data/fish6/asm/verkko_asm_v1.4.1/assembly.fasta fish6_index

#bowtie-build /data/okendojo/zebrafish/data/g3/assembly/verkko_asm_run2/assembly.fasta g3_asm_index
