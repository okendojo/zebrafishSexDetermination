#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=210g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=F11_sarek_variantcall

#load modules
module add nextflow
module add singularity


cd  /data/okendojo/zebrafish/data/g3/F2_variants/Fish11

fasta="/data/okendojo/zebrafish/data/fish11/finalPolish/nhgri_polished_fish11.fasta"
fai="/data/okendojo/zebrafish/data/fish11/finalPolish/nhgri_polished_fish11.fasta.fai"

nextflow run nf-core/sarek --input 'samplesheet.csv'  --genome null  --igenomes_ignore  --fasta $fasta   --fasta_fai $fai  --outdir 'nhgri_fish11'   --tools 'manta'  --trim_fastq true   --concatenate_vcfs false  --trim_fastq true  --skip_tools baserecalibrator  -resume  -profile biowulf
