#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=fastqTofasta

#load modules

module load seqtk
cd /data/okendojo/zebrafish/data/g3/pacBio

seqtk seq -A /data/okendojo/zebrafish/data/g3/pacBio/hifi.fastq.gz > hifi.fasta 
