#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=210g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=G3_correction

#load modules
module load herro

cd /data/okendojo/zebrafish/data/AB 

herro inference -m model_R10_v0.1.pt -b 64 -t 16 /data/okendojo/zebrafish/data/g3/ontData/ONT.fastq.gz /data/okendojo/zebrafish/data/g3/g3_ont_herro.fasta
