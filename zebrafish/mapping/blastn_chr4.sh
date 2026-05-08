#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=f6blast

#Load the module
module load blast

cd /data/okendojo/zebrafish/data/AB/polishing/rDNAcopynumber

#makeblastdb -in ../asm_verkko2_2/assembly.fasta -dbtype nucl -out db

#blastn -query /data/okendojo/zebrafish/data/AB/batches_ont/ont_ul_100kfilt.fasta -db rdna.fa -out AB_rDNA_hits.tsv -outfmt 6
blastn -query  /data/okendojo/zebrafish/data/fish6/ontData/concatenated/ont_100kb.fasta -db rdna.fa -out fish6_rDNA_hits2.tsv -outfmt 6
#blastn -query /data/okendojo/zebrafish/data/fish6/asm/tangleResolution/rDNA.fasta -db db -out results.txt -outfmt 6

