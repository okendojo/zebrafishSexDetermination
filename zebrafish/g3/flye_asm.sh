#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=90:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=fly_asm

#load modules
module add flye

cd /data/okendojo/zebrafish/data/g3

#GIZ
flye --nano-hq /data/okendojo/zebrafish/data/g3/ontData/individual/*.fastq.gz  --out-dir ont_asm --threads 32

