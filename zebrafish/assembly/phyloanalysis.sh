#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=phyloanalysis


#==================================
#####Load the required module======
#==================================
module load R


#move to the right dir

cd /data/okendojo/zebrafish/data/fish6/centromere 

Rscript phylo.R 
