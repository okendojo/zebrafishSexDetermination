#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=mkref-rna

#Load the module
module load cellranger
module load cellranger-arc

#Download the zebrafish ref from: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_mr


cd /data/okendojo/scRNA_project/hui_data/rawData/fastqs_GEX


cellranger mkref \
  --genome=GRCz11 \
  --fasta=/data/okendojo/zebrafish/refGenome/ref/Danio_rerio.GRCz11.dna.primary_assembly.fa \
  --genes=/data/okendojo/zebrafish/refGenome/ref/Danio_rerio.GRCz11.105_filtered.gtf --nthreads=24 --memgb=200
