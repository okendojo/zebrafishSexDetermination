#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=90:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=splice_quant

#load modules
module add nextflow
module add singularity


cd /data/okendojo/zebrafish/data/g3/eQTL/splicequant

nextflow run nf-core/rnasplice --input 'rnavar_samplesheet.csv' --outdir 'splice_result' -profile biowulf  --contrasts "contrastsheet.csv" -resume --fasta '/data/okendojo/zebrafish/data/g3/rna_sequences/gene_quantification/ncbi_dataset/Danio_rerio.GRCz11.dna.primary_assembly.fa' --gtf '/data/okendojo/zebrafish/data/g3/rna_sequences/gene_quantification/ncbi_dataset/Danio_rerio.GRCz11.110.chr.gtf' --star_index '/data/okendojo/zebrafish/refGenome/star_index'   --skip_bigwig false   --skip_fastqc true  --rmats true  --rmats_novel_splice_site true  --dexseq_exon true  --edger_exon true  --dexseq_dtu false  --suppa false  --sashimi_plot true   --aligner star_salmon  
