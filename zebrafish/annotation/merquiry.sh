#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=80g
#SBATCH --ntasks-per-core=1
#SBATCH --time=96:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=merquiry

#move to the working directory

#Load the modules
module load python/3.9
module load snakemake/7.7.0
module load R/4.2.0
module load bedtools/2.30.0
module load samtools/1.9

cd /data/okendojo/zebrafish/data/g3/rna_sequences/blob_plot

$tools/merqury/merqury.sh /data/okendojo/zebrafish/data/g3/illumina/meryldbs/child_RNA.meryl /data/okendojo/zebrafish/data/g3/illumina/meryldbs/WIK.meryl /data/okendojo/zebrafish/data/g3/illumina/meryldbs/TL.meryl /data/okendojo/zebrafish/data/g3/rna_sequences/trimmed/transcript.fasta wik_tl
