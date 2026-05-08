#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=210g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=isoseq

#load modules
module add nextflow
module add singularity

cd /data/okendojo/zebrafish/data/annotation/fish6/annot

nextflow run nf-core/isoseq -profile biowulf --input "/data/okendojo/zebrafish/data/annotation/fish6/samplesheet.csv" --outdir isoseqannot --fasta "/data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/asm_polished/NHGRI_Fish6_cons.fasta" --primers primer.fasta --gtf "/data/okendojo/zebrafish/data/fish6/liftoff/fish6_asm/fish6_annotation.gtf" --aligner minimap2 	
