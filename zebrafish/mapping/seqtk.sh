#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=20:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=seqtk

#Load the module

ml seqtk

cd /data/okendojo/zebrafish/data/AB/polishing/rDNAcopynumber

seqtk subseq /data/okendojo/zebrafish/data/AB/batches_ont/ont_ul.fasta reads.txt > abrDNA.fasta
