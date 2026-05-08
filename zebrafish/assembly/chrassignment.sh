#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=chr_assignment


#Load the modules
source /data/$USER/conda/etc/profile.d/conda.sh && source /data/$USER/conda/etc/profile.d/mamba.sh

#activate ragtag
mamba activate ragtag

cd /data/okendojo/zebrafish/data/fish6/asm/asm_verkko_24x

ragtag.py scaffold GRCz11.fasta assembly.fasta -w -u  -t 24 -o chromosome_assignment/

