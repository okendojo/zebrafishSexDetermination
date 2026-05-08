#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=sckind

#load modules
module load meryl/1.4.1 snakemake/7.32.4 pigz R samtools


cd /data/okendojo/zebrafish/data/g3/sex_project/svbyeye/scikind/1ZZ_chr4

snakemake -c 24 -s SCINKD/SCINKD.v2.1.0.FULL.snakefile

