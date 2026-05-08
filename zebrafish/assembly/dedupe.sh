#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --gres=lscratch:500
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=dedupe

export TMPDIR=/lscratch/$SLURM_JOB_ID

cd /data/okendojo/zebrafish/data/fish11/ont/fastq_pass/concatenatedONT

gunzip zfish11ONT.fastq.gz

czid-dedup -i zfish11ONT.fastq -o dedupedONT.fastq -c custom-cluster.csv 
