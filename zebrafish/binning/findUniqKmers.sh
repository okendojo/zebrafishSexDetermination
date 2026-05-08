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
module load python
module add gcc kmc
module load openmpi/4.1.3/gcc-11.3.0



cd /data/okendojo/zebrafish/data/ab_asm/binning
mkdir findUniqKmers

cd findUniqKmers

find-unique-kmers -k 21 -p 24 /data/okendojo/zebrafish/data/g3/illumina/AB_F_CB/HFVYJDSX5_19412535_S54_L002_R1_001.fastq.gz,/data/okendojo/zebrafish/data/g3/illumina/AB_F_CB/HFVYJDSX5_19412535_S54_L002_R2_001.fastq.gz \
       /data/okendojo/zebrafish/data/g3/illumina/TU_M_CB/HFVYJDSX5_19412541_S50_L002_R1_001.fastq.gz,/data/okendojo/zebrafish/data/g3/illumina/TU_M_CB/HFVYJDSX5_19412541_S50_L002_R2_001.fastq.gz 	
