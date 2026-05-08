#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=addtelo_scaffold

#load modules
module load ragtag
module load seqkit

cd /data/okendojo/zebrafish/data/fish11/polishing/telomere

python addTelomere.py
