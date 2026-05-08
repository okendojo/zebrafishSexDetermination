#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=90:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=findUniqKmers

#activate conda
conda activate trio_binning

#load modules
#module load openmpi/4.1.3/gcc-11.3.0
#load modules
#module add gcc #kmc

cd /data/okendojo/zebrafish/data/g3/binning/TU_AB

find-unique-kmers -k 21 -p 24 HFVYJDSX5_19412535_S54_L002_R1_001_trimmed.fq.gz,HFVYJDSX5_19412535_S54_L002_R2_001_trimmed.fq.gz HFVYJDSX5_19412541_S50_L002_R1_001_trimmed.fq.gz,HFVYJDSX5_19412541_S50_L002_R2_001_trimmed.fq.gz


cd /data/okendojo/zebrafish/data/g3/binning/TL_WIK

find-unique-kmers -k 21 -p 24 HFVYJDSX5_19412527_S51_L002_R1_001_trimmed.fq.gz,HFVYJDSX5_19412527_S51_L002_R2_001_trimmed.fq.gz HFVYJDSX5_19412539_S55_L002_R1_001_trimmed.fq.gz,HFVYJDSX5_19412539_S55_L002_R2_001_trimmed.fq.gz 


