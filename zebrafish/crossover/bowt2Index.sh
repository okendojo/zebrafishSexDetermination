#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=bowtie2_index

cd /data/okendojo/zebrafish/refGenome/index

#Load the modules
module load bowtie/2-2.3.4.3

bowtie2-build --threads 32 GRCz11_genomic.fasta index
