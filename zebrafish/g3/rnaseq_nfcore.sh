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
module load nextflow
module load singularity


cd /data/okendojo/zebrafish/data/g3/sex_project/asm/t2tMap/bam_fixed_rg/na_rna_fastq

nextflow run nf-core/rnaseq --input '/data/okendojo/zebrafish/data/g3/sex_project/asm/t2tMap/bam_fixed_rg/na_rna_fastq/fastq/samplesheet.csv' --outdir 'nadia_rna_quant' -profile biowulf  --email 'javanokendo@gmail.com' --with_umi false -resume --fasta '/data/okendojo/zebrafish/data/g3/sex_project/asm/t2tMap/GCF_049306965.1_GRCz12tu_genomic.fna' --gtf '/data/okendojo/zebrafish/data/g3/sex_project/asm/t2tMap/GCF_049306965.1_GRCz12tu_genomic.gtf'   --skip_qc false  --skip_stringtie false --skip_fastqc false --skip_preseq false --skip_dupradar false --skip_rseqc false --skip_biotype_qc false --skip_deseq2_qc false --skip_multiqc false --skip_bigwig false   --with_umi false   --skip_trimming false 
