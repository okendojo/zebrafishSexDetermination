#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=mafft


#Load the modules
module load mafft


cd /data/okendojo/zebrafish/data/fish6/centromere/anal_phylo

mafft  f11.fasta > f11_aligned.fasta 
