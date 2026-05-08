#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=Select_stainGlass

cd /data/okendojo/zebrafish/data/AB/asm/polishing/stainglass/StainedGlass

#load modules
module load python
#source /data/$USER/conda/etc/profile.d/conda.sh && source /data/$USER/conda/etc/profile.d/mamba.sh
#mamba activate snakemake

module load snakemake/8.16.0
module load minimap2
module load samtools
module load bedtools
#module load R

time snakemake --cores 32 --config sample=NOTAL fasta=/data/okendojo/zebrafish/data/AB/asm/polishing/stainglass/largeUnmappedTU11Tu12.fasta

snakemake --cores 24 make_figures
