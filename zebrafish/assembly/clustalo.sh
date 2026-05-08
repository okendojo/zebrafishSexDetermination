#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=clustalo


#Load the modules
module load clustalo

cd /data/okendojo/zebrafish/data/fish6/centromere/anal_phylo

clustalo -i f11.fasta -o f11_6_alig.fasta --threads=24 -v
