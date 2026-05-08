#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=lastz


#Load the modules
module load lastz

cd /data/okendojo/zebrafish/data/fish6/centromere

lastz fish11_centromera.fasta[multiple] --self --format=general  --output=f11_alignment.txt --chain  --gapped  --inner=2000 --ydrop=3400
