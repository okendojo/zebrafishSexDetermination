#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=210g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=cage_cnvs

#load modules
module add nextflow
module add singularity


cd  /data/okendojo/zebrafish/data/g3/mir_430Gene/cage_data


nextflow run nf-core/sarek --input '/data/okendojo/zebrafish/data/g3/mir_430Gene/cage_data/samplesheet/danio_samplesheet.csv'  --genome null  --igenomes_ignore  --fasta '/data/okendojo/zebrafish/data/g3/mir_430Gene/cage_data/assembly/chr4.fasta'   --fasta_fai '/data/okendojo/zebrafish/data/g3/mir_430Gene/cage_data/assembly/chr4.fasta.fai' --outdir 'cage_cvns' --tools 'cnvkit'   --save_output_as_bam true --trim_fastq true  --skip_tools baserecalibrator  -resume  -profile biowulf
