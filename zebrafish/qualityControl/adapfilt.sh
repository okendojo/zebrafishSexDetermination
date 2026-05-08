#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=F6_Filt


#==================================
#####Load the required module======
#==================================
module load bamtools
module load blast

cd /data/okendojo/paradisfishProject/child

mkdir qc

hifiadapterfilt.sh -o qc 
