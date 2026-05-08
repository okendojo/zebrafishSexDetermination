#!/bin/sh
################
## HEADER FOR BIOWULF BASH SCRIPT
################
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=750g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=G3_hifiasm_ont
###
## MODULES TO USE
###
module load hifiasm
module load seqkit
##########
## BEGIN SCRIPT
##########
export TMPDIR=/lscratch/$SLURM_JOB_ID

cd  /data/okendojo/zebrafish/data/g3/assembly/hifiasm

hifiasm -t 24 -o ont_hifi_asm --ul /data/okendojo/zebrafish/data/g3/ontData/ont_uniq.fastq.gz /data/okendojo/zebrafish/data/g3/pacBio/*.fastq.gz
