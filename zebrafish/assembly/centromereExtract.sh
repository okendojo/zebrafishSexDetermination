#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=centrExtract


#Load the modules
ml jellyfish tandem-genotypes samtools bedtools seqtk

cd  /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/mapping/ont

#samtools view --threads 24 -O BAM  --min-MQ 20 -b -L chr17cent.bed  ONT.pri.bam -o extracted_cent.bam

#bedtools bamtofastq -i extracted_cent.bam -fq extracted_cent.fastq
#seqtk seq -A extracted_cent.fastq > extracted_cent.fasta

#tandemmapper.py --nano extracted_cent.fasta chr17.fasta -t 24 -o chrom17centres

tandemquast.py  --nano extracted_cent.fasta chr17.fasta -t 24 -o chrom17centres
