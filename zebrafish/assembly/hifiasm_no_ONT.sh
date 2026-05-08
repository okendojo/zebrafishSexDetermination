#!/bin/sh
################
## HEADER FOR BIOWULF BASH SCRIPT
################
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=750g
#SBATCH --ntasks-per-core=3
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=hifiasm_ont
###
## MODULES TO USE
###
module load hifiasm
##########
## BEGIN SCRIPT
##########
export TMPDIR=/lscratch/$SLURM_JOB_ID

cd  /data/okendojo/zebrafish/data/AB/asm/hifiasm_ONT

hifiasm -t 30 --ont  -o ab_asm /data/okendojo/zebrafish/data/AB/batches_ont/ont_ul.fastq.gz

#hifiasm -t 30 -o ont_hifi_asm --ul /data/okendojo/zebrafish/data/AB/batches_ont/ont_ul.fastq.gz /data/okendojo/zebrafish/data/AB/hifi/*.fastq.gz
