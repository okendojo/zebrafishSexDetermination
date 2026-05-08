#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --gres=lscratch:500
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=trim_f6


#==================================
#####Load the required module======
#==================================
export TMPDIR=/lscratch/$SLURM_JOB_ID

module load bamtools
module load porechop
module load blast

cd /data/okendojo/zebrafish/data/fish6/ontData/untrimmed

mkdir ztrimmed

for file in *.fastq.gz

do

porechop -i "${file}" -o ztrimmed/"${file%.fq}.trim.fastq.gz" -t 32

done
