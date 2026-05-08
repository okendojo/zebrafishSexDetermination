#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=360g
#SBATCH --ntasks-per-core=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=cent_staineGlass

#load module
#module load python

#source /data/$USER/conda/etc/profile.d/conda.sh 
#module load mamba_install
source myconda
#mamba_install --init-file=~/bin/myconda --init-only
conda activate snakemake

cd /data/okendojo/zebrafish/data/AB/polishing/centrominer/stainedglass/StainedGlass

snakemake --cores 24  make_figures

