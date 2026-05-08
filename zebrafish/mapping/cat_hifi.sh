#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=cat_hifi

cd /data/okendojo/zebrafish/data/fish11/hifi

cat m64467e_230509_164640.hifi_reads.fastq.gz  m64467e_230511_034148.hifi_reads.fastq.gz > hifi.fastq.gz

gunzip  hifi.fastq.gz
