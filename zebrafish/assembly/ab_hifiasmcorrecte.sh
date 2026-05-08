#!/bin/sh
################
## HEADER FOR BIOWULF BASH SCRIPT
################
#SBATCH --job-name=F6_T2T
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 48
#SBATCH --partition=largemem
#SBATCH --mem=800G
#SBATCH --gres=lscratch:500
#SBATCH --time=240:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov

###
## MODULES TO USE
###
module load hifiasm
##########
## BEGIN SCRIPT
##########
export TMPDIR=/lscratch/$SLURM_JOB_ID

cd /data/okendojo/zebrafish/data/AB/asm/hifiasm

hifiasm -o AB_hifiasm_ul.asm -t 48 --telo-m CCCTAA -l 0 /data/okendojo/zebrafish/data/AB/hifi/*.fastq.gz --ul /data/okendojo/zebrafish/data/AB/prelimONT/*.fastq.gz

#hifiasm -o AB_hifiasm_ul.asm -t 48 -l 0 --ul `ls *.fastq.gz | tr '\n' ',' | perl -ple 'chop'` /data/okendojo/zebrafish/data/fish6/hifi/*.fastq.gz

