#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=combineONT


#Load the modules
module load seqtk

cd /data/okendojo/zebrafish/data/fish11/ont

cat ONT/*.fastq.gz > concatenate/ont_fastq.gz

seqtk seq -A concatenate/ont_fastq.gz > concatenate/ont.fasta

seqtk hpc concatenate/ont.fasta > concatenate/ont_hpc.fasta

grep ">" concatenate/ont.fasta | sort | uniq -c | sed 's/>//' | awk -v OFS="\t" '{ print $1, $2}' | awk '$1>1' > concatenate/ont_dup_entries.txt
