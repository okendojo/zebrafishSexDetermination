#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=512g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=G3_hifiasm

#Hifiasm without binning with 
module load hifiasm

#Normal assembly without triobinning;# Assemble heterozygous genomes with built-in duplication purging

cd /data/okendojo/zebrafish/data/g3/assembly

hifiasm -o hifiasm.asm -t 32  -l 1  /data/okendojo/zebrafish/data/g3/pacBio/*.fastq.gz --ul /data/okendojo/zebrafish/data/g3/ontData/combined/ont.fastq.gz


