#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=seqtk


#Load the modules
module load seqtk
module load seqkit


cd /data/okendojo/zebrafish/data/g3/ontData

#cat *.fastq.gz | seqtk seq -L 100000 - | gzip - > zont_100kb.fastq.gz

#seqtk subseq ont.fasta ont_subset.id > ont_subset.fasta

#cat /data/okendojo/zebrafish/data/g3/ontData/individual/*.fastq.gz | seqtk seq -A  - > comb_ont.fasta 

#grep ">" comb_ont.fasta | awk -v OFS="\t" '{ print $1}' | sed 's/>//g' | sort | uniq -c | awk -v OFS="\t" '{ print $1, $2}' > comb.txt 

#gzip -k ont_subset.fasta 

#seqtk seq -A ont.fastq.gz | seqkit rmdup -s -j 20 - | seqkit seq -m 10000 - > ont2.fasta 

#seqkit seq -m 10000 ont.fastq.gz | seqkit rmdup -s -j 20 -  > ont_filt10k.fastq.gz

seqkit rmdup ont.fastq.gz -o ont_uniq.fastq.gz -D tmpdup.txt -d tmpdup.fastq.gz 
