#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=360g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=ABG3


#Load the modules
module load verkko
module load snakemake

#Run the reads extraction first

cd  /data/okendojo/zebrafish/data/AB/asm

#Run verkko
verkko -d AB_G3AB --hifi /data/okendojo/zebrafish/data/AB/hifi/*.fastq.gz --nano /data/okendojo/zebrafish/data/AB/AB_G3AB.fasta --slurm

#verkko -d AB_G3AB_tangleresolved --assembly AB_G3AB --cns-run 10 90 90  --hifi /data/okendojo/zebrafish/data/AB/hifi/*.fastq.gz --nano /data/okendojo/zebrafish/data/AB/AB_G3AB.fasta --slurm --paths /data/okendojo/zebrafish/data/AB/asm/polishing/tangleresolution/paths.gaf
