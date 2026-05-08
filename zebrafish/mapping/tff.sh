#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=trf

#load the modules

module load minimap2
module load syri
module load plotsr
module load trf

cd /data/okendojo/zebrafish/data/fish6/asm/variationAnalysis

trf fish11.fasta 2 7 7 80 10 50 500 -f -d -m

trf fish6.fasta 2 7 7 80 10 50 500 -f -d -m
