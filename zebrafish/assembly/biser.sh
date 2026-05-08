#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=biser2


cd /data/okendojo/zebrafish/data/AB/polishing/segmentaldups

#Run repeatmasker
ml repeatmasker


# Inputs
GENOME="/data/okendojo/zebrafish/data/AB/asm/ab_genome_assembly/GRCz12ab.fasta"                # your assembly FASTA
OUTDIR="RM_out"
THREADS=16

mkdir -p "$OUTDIR"

RepeatMasker -e ncbi -pa $(expr $SLURM_CPUS_PER_TASK / 2)  -xsmall -gff -dir "$OUTDIR" -species "vertebrates" ${GENOME}


#load modules
module load samtools

biser -t 30 --max-error 20 --max-edit-error 10 --output grcz12ab_sd.sd --gc-heap 3G --kmer-size 31 RM_out/GRCz12ab.fasta.masked 
