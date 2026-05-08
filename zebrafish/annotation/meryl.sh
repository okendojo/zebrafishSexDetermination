#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=170:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=maryl

#Load the required modules
module add maryl

cd /data/okendojo/zebrafish/data/g3/rna_sequences/trimmed

meryl count k=21 threads=24 *.fq.gz output child_RNA.meryl
