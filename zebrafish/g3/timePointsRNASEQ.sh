#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=rnaseq_quant

#load modules
module add nextflow
module add singularity


cd /data/okendojo/zebrafish/data/g3/rna_sequences/timepointsQuant

nextflow pull nf-core/rnaseq

nextflow run nf-core/rnaseq --input '/data/okendojo/zebrafish/data/g3/rna_sequences/untrimmed/samplesheet_2.csv' --outdir 'binned_quant' -profile biowulf --email 'javanokendo@gmail.com' --with_umi false -resume --fasta '/data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.dna.primary_assembly.fa' --gtf '/data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.111.gtf'   --remove_ribo_rna true  --skip_qc false  --skip_stringtie false --skip_fastqc false --skip_preseq false --skip_dupradar false --skip_rseqc false --skip_biotype_qc false --skip_deseq2_qc false --skip_multiqc false --skip_bigwig false   --with_umi false   --skip_trimming false 
