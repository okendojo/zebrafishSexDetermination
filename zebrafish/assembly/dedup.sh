#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=fastp


#Load the modules
module load fastp

cd /data/okendojo/zebrafish/data/fish6/ontData/concatenated 

fastp -i  f6_ont.fastq.gz -o f6_ont_deduped.fastq.gz --dedup --length_required 10000 --disable_adapter_trimming --compression 4 --thread 24  
