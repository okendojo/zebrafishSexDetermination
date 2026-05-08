#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=1000g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=seqtk

# Load meryl module
module add seqtk

cd /data/okendojo/zebrafish/data/g3/rna_sequences/haplotype_gene_quantification

seqtk subseq forward_reads_R1.fastq.gz txt/AB_R1.txt > AB_R1.fastq 

seqtk subseq forward_reads_R1.fastq.gz txt/TU_R1.txt > TU_R1.fastq

seqtk subseq forward_reads_R1.fastq.gz txt/TL_R1.txt > TL_R1.fastq

seqtk subseq forward_reads_R1.fastq.gz txt/WIK_R1.txt > WIK_R1.fastq


#seqtk subseq reverse_reads_R2.fastq.gz AB_R2.txt > AB_R2.fastq 

#seqtk subseq reverse_reads_R2.fastq.gz TU_R2.txt > TU_R2.fastq

#seqtk subseq reverse_reads_R2.fastq.gz TL_R2.txt > TL_R2.fastq

#seqtk subseq reverse_reads_R2.fastq.gz WIK_R2.txt > WIK_R2.fastq
