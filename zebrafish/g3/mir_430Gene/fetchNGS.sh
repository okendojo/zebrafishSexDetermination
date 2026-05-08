#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=90:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=fetchCAGE

#load modules
module add nextflow
module add singularity

cd /data/okendojo/zebrafish/data/g3/mir_430Gene

nextflow run nf-core/fetchngs -profile biowulf --input cage_accession.csv --outdir cage_data --download_method sratools 
