#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=verkko_consensus

#load module 
module add python/3.9

cd /data/okendojo/zebrafish/data/fish11/gapfilling

python3 gaps2bed.py GRCz11.fasta grcz_gaps.bed
