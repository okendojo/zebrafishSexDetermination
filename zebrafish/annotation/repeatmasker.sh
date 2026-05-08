#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=170:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=ZFish_repmasker

#Load the required modules
module add repeatmasker


cd /data/okendojo/zebrafish/data/fish11/annotation

# Define paths and files
genome_fasta="fish11_merfin_polished.fasta"         # Path to your genome assembly in FASTA format


# Run RepeatMasker to mask repeats in the genome assembly
RepeatMasker -pa 24 -species "danio rerio" -xsmall -html -s -gff -e ncbi -dir repResults "$genome_fasta"
