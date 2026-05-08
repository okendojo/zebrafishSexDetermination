#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=230g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=fastp


ml fastp

cd  /data/okendojo/zebrafish/data/g3/illumina/TU_M_CB

fastp -i HFVYJDSX5_19412541_S50_L002_R1_001.fastq.gz -I HFVYJDSX5_19412541_S50_L002_R2_001.fastq.gz -o merged_R1.fastq.gz -O merged_R2.fastq.gz --merged_out TU.fastq.gz --merge

cd /data/okendojo/zebrafish/data/g3/illumina/AB_F_CB

fastp -i HFVYJDSX5_19412535_S54_L002_R1_001.fastq.gz -I HFVYJDSX5_19412535_S54_L002_R2_001.fastq.gz -o merged_R1.fastq.gz -O merged_R2.fastq.gz --merged_out AB.fastq.gz --merge

cd /data/okendojo/zebrafish/data/g3/illumina/TL_F_CB

fastp -i HFVYJDSX5_19412527_S51_L002_R1_001.fastq.gz  -I HFVYJDSX5_19412527_S51_L002_R2_001.fastq.gz -o merged_R1.fastq.gz -O merged_R2.fastq.gz --merged_out TL.fastq.gz --merge

cd /data/okendojo/zebrafish/data/g3/illumina/WIK_M_CB

fastp -i HFVYJDSX5_19412539_S55_L002_R1_001.fastq.gz -I HFVYJDSX5_19412539_S55_L002_R2_001.fastq.gz -o merged_R1.fastq.gz -O merged_R2.fastq.gz --merged_out WIK.fastq.gz --merge


