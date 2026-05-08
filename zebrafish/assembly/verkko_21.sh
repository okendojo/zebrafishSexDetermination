#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=CBZW_9_asm_verkko


#Load the modules
ml verkko
ml python/3.10

#move to the wkdir
cd  /data/okendojo/zebrafish/data/g3/sex_project

#Run verkko
verkko -d CBZW_9_asm --hifi /data/okendojo/zebrafish/data/g3/sex_project/CBZW_9/m84270_250710_154513_s2.hifi_reads.B11_bc2082.fastq --no-nano --slurm --lay-run 16 100 96
