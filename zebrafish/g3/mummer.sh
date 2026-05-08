#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=nucmer_run

#Load the modules
module load mummer

cd /data/okendojo/zebrafish/data/g3/assembly/verkko_asm2

nucmer --mum -t 24 -p AB_haplotype -b 200 -c  100 GRCz11.fasta AB.fasta 

nucmer --mum -t 24 -p TU_haplotype -b 200 -c 100 GRCz11.fasta TU.fasta

nucmer --mum -t 24 -p TL_haplotype -b 200 -c 100 GRCz11.fasta TL.fasta

nucmer --mum -t 24 -p WIK_haplotype -b 200 -c 100 GRCz11.fasta WIK.fasta 
