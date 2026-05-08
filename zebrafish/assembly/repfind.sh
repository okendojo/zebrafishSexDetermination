#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=seqtk


#Load the modules
module load seqkit
module load R

cd /data/okendojo/zebrafish/data/fish6/centromere

sh Setup_Run_Repeats.sh -i zebrafish -f fish11.fasta -c 30  -g TRUE
