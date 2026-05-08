#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=750g
#SBATCH --ntasks-per-core=1
#SBATCH --time=140:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=nadia_asm

#load modules
ml hifiasm

cd /data/okendojo/zebrafish/data/g3/sex_project/asm

hifiasm -o nadia_asm -t 24 /data/okendojo/zebrafish/data/g3/sex_project/fastq_files/nadia/*.fastq.gz


