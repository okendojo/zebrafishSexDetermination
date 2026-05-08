#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=fish11_2_t2t


#Load the modules
module load verkko/2.2
module load seqtk


#Run the reads extraction first
cd /data/okendojo/zebrafish/data/fish11

#Run verkko
verkko --assembly asm_verkko2_2 -d chrs_t2t  --hifi /data/okendojo/zebrafish/data/fish11/hifi/*.fastq.gz --nano /data/okendojo/zebrafish/data/fish11/ont/concatenate/ont_fastq.gz  --paths /data/okendojo/zebrafish/data/fish11/polishing/consensus/tmp.gaf --slurm 
