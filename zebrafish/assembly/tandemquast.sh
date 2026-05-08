#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=tandemQuast


#Load the modules
module load jellyfish

cd /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/nucfreq/chromosomes


for file in *.fasta ; do tandemquast.py -t 24 -o ${file}_quast --nano /data/okendojo/zebrafish/data/fish6/ontData/concatenated/ont_filt.fasta ${file} ; echo "Now working on $fasta" ; done 
