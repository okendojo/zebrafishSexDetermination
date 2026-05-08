#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=750g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --gres=lscratch:100
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=HOR22

#load module
ml R trf

cd /data/okendojo/zebrafish/data/AB/polishing/centrominer/splitcent

Rscript hor_from_cent.R /data/okendojo/zebrafish/data/AB/polishing/syri/AB.fa myHORs

#trf /data/okendojo/zebrafish/data/AB/polishing/syri/AB.fa 2 5 7 80 10 50 500 -d -f
