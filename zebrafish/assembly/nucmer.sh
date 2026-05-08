#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=10:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=mashmap


#==================================
#####Load the required module======
#==================================
module load mummer

cd /data/okendojo/zebrafish/data/segmental_duplication/synteny/fasta

nucmer --prefix=asm_comparison F6.fasta F11.fasta 
delta-filter -r -q asm_comparison.delta > asm_comparison.filtered.delta
show-coords -THrd asm_comparison.filtered.delta > structural_variants.txt

