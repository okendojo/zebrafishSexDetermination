#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=500g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=var_exon_extraction


#Load the modules
module load python

cd /data/okendojo/zebrafish/data/g3/F1_variants/phasing_variants

python ana2.py
