#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=800g
#SBATCH --gres=lscratch:500
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=porechopONT


#==================================
#####Load the required module======
#==================================
export TMPDIR=/lscratch/$SLURM_JOB_ID

module load bamtools
module load porechop
module load blast

cd /data/okendojo/zebrafish/data/fish11/ont/fastq_pass/zfishcomb

porechop -i zfishcombined.fastq.gz -o fish11ontTrimmed.fastq.gz --threads 32


