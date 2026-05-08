#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=90:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=fish11_quant

#load modules
module add nextflow
module add singularity


cd /data/okendojo/zebrafish/data/g3/mir_430Gene/cage_data/fish11_mirquant

nextflow run nf-core/rnaseq --input '/data/okendojo/zebrafish/data/g3/mir_430Gene/cage_data/samplesheet/danio_samplesheet.csv' --outdir 'fish11_mir430'  -profile biowulf --email 'javanokendo@gmail.com' --with_umi false -resume --fasta '/data/okendojo/zebrafish/data/g3/mir_430Gene/cage_data/assembly/fish11.fa' --gtf '/data/okendojo/zebrafish/data/g3/mir_430Gene/cage_data/assembly/genes.gtf'  --save_unaligned false --skip_stringtie false --skip_fastqc false --skip_preseq false --skip_dupradar false --skip_rseqc false --skip_biotype_qc false --skip_deseq2_qc false --skip_multiqc false --skip_bigwig false
