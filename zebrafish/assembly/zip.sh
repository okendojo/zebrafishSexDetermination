#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=zip


#Load the modules
module load seqtk

cd /data/okendojo/zebrafish/data/AB

gzip ont_filtered.fastq

seqtk seq -A ont_filtered.fastq.gz > ont_filtered.fasta
