#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=751g
#SBATCH --ntasks-per-core=1
#SBATCH --time=10:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=100kasm_seqkit


#Load the modules
module load seqkit
module load seqtk
module load hifiasm

cd /data/okendojo/paradisfishProject/child/ont

seqtk seq -A mop.nanopore.pass.fastq.gz > ont.fasta

seqkit rmdup -s -j 24  ont.fasta -o ont_deduped.fasta 
