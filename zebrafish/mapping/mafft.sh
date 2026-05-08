#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=mafft

#load the modules
module load mafft

cd /data/okendojo/zebrafish/data/centromere_analysis
 
mafft --auto centromere_sequences.fasta > f11_mafft_output.aln

