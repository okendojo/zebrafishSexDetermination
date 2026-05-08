#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=cellRanger

#Load the required modules
module purge
module load nextflow
module load singularity


#move to the directory containing the data
cd /data/okendojo/scRNA_project/hui_data


CONFIG=/data/okendojo/slugProject/annot/scrna/nextflow.config
FASTA='/data/okendojo/zebrafish/refGenome/GRCz11_genomic.fasta'
GTF='/data/okendojo/zebrafish/refGenome/GCF_000002035.6_GRCz11_genomic.gtf'

mkdir results

nextflow run nf-core/scrnaseq --input samplesheet.csv -c ${CONFIG} --outdir results -resume --fasta '/data/okendojo/zebrafish/refGenome/GRCz11_genomic.fasta' --gtf '/data/okendojo/zebrafish/refGenome/GCF_000002035.6_GRCz11_genomic.gtf' --protocol 10XV2 --aligner star -profile singularity
