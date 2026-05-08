#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=200g
#SBATCH --ntasks-per-core=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=seqkit


#Load the modules
module load seqkit

cd /data/okendojo/zebrafish/data/fish6/ontData/concatenated

seqkit split2 -p 10 f6_ont.fastq.gz --threads 24

