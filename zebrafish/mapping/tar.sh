#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=TAR_FILE

cd /data/okendojo/zebrafish/data/fish11/ont/fastq_pass/concatenatedONT

gzip zfish11ONT.fastq

tar -czvf ont_fastq.tar.gz zfish11ONT.fastq.gz


