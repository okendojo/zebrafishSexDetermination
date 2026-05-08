#!/bin/bash 
#BATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=bwa_index

cd /data/okendojo/zebrafish/refGenome

module load bwa 

bwa index GRCz11_genomic.fasta
