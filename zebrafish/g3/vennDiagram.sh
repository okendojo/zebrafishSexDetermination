#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=1000g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=vennDiagram


#Load the modules
module add R

#move to the dir containing the data

cd /data/okendojo/zebrafish/data/g3/rna_sequences/meryl_analysis

Rscript vennD.R
