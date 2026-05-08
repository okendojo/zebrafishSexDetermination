#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=210g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=F6_sarek_variantcall

#load modules
module add nextflow
module add singularity


cd  /data/okendojo/zebrafish/data/segmental_duplication/synteny/SV

fasta="/data/okendojo/zebrafish/data/segmental_duplication/synteny/fasta/F6.fasta"
fai="/data/okendojo/zebrafish/data/segmental_duplication/synteny/fasta/F6.fasta.fai"

nextflow run nf-core/sarek --input 'samplesheet.csv'  --genome null  --igenomes_ignore  --step variant_calling --fasta $fasta   --fasta_fai $fai  --outdir 'nhgrifish11_fish6'  --tools 'manta'  --trim_fastq true   --concatenate_vcfs false  --trim_fastq true  --skip_tools baserecalibrator  -resume  -profile biowulf
