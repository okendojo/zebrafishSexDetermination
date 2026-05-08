#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=170:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=centromere_analysis


#Load the required modules
module add k8

DIR="/data/okendojo/zebrafish/data/centromere_analysis"

cd ${DIR}

FISH11="/data/okendojo/zebrafish/data/centromere_analysis/fish11.fasta"
FISH6="/data/okendojo/zebrafish/data/centromere_analysis/fish6.fasta"
GRCZ11="/data/okendojo/zebrafish/data/centromere_analysis/GRCz11.fasta"

#=====================================================
#===========START FISH11 ANALYIS======================
#=====================================================

# Run RepeatMasker to mask repeats in the genome assembly
RepeatMasker -pa 24 -species "danio" -xsmall -small -html -s -gff -e ncbi -dir fish11_repeats "$FISH11"

# convert .out file to bedfile
k8 dna-nn/parse-rm.js fish11_repeats/fish11.fasta.out > fish11_rm.bed

# Generate training data in FASTQ. Base qualities indicate labels.
gen-fq -m2 ${FISH11} fish11_rm.bed > fish11.fq

# Training. We trained 10 models with different random seeds
dna-brnn -t8 -n32 -b5m -m50 -s14 -o fish11-alpha.knm fish11.fq

# find (ATTCC)n and alpha satellites for long contigs
## dna-brnn -Ai fish11-alpha.knm -t16 ${FISH11} > fish11_seq.bed

#=======================================================
#=============START FISH6 ANALYSIS======================
#=======================================================
RepeatMasker -pa 24 -species "danio" -xsmall -small -html -s -gff -e ncbi -dir fish6_repeats "$FISH6"

k8 dna-nn/parse-rm.js fish6_repeats/fish6.fasta.out > fish6_rm.bed

gen-fq -m2 ${FISH6} fish6_rm.bed > fish6.fq

dna-brnn -t8 -n32 -b5m -m50 -s14 -o fish6-alpha.knm fish6.fq

#Train the model
## dna-brnn -Ai fish6-alpha.knm -t16 ${FISH6} > fish6_seq.bed


#Run repeatmasker using GRCz11
RepeatMasker -pa 24 -species "danio" -xsmall -small -html -s -gff -e ncbi -dir GRCZ11_repeats "$GRCZ11"


