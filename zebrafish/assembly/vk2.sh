#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=AB_vk2


#Load the modules
module load verkko

#Run the reads extraction first
cd  /data/okendojo/zebrafish/data/AB/asm

#Run verkko
verkko -d  asm_comp250505 --hifi /data/okendojo/zebrafish/data/AB/hifi/*.fastq.gz --nano /data/okendojo/zebrafish/data/AB/ont_files/*.fastq.gz --slurm  --lay-run 16 100 96 --spl-run 16 100 96 --mbg-run 16 100 96
