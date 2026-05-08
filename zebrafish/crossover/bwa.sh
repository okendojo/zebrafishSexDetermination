#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=bwa_ZF

cd /data/okendojo/zebrafish/data/G3/TU_AB_alignment

#Load the modules
module load bwa
module load samtools

#See the data paths
REF_GENOME=/data/okendojo/zebrafish/refGenome/GRCz11_genomic.fasta
TU_M_R1=/data/okendojo/zebrafish/data/G3/TU_M_CB/qc_results/HFVYJDSX5_19412541_S50_L002_R1_001_trimmed.fq.gz
TU_M_R2=/data/okendojo/zebrafish/data/G3/TU_M_CB/qc_results/HFVYJDSX5_19412541_S50_L002_R2_001_trimmed.fq.gz
AB_F_R1=/data/okendojo/zebrafish/data/G3/AB_F_CB/qc_results/HFVYJDSX5_19412535_S54_L002_R1_001_trimmed.fq.gz
AB_F_R2=/data/okendojo/zebrafish/data/G3/AB_F_CB/qc_results/HFVYJDSX5_19412535_S54_L002_R2_001_trimmed.fq.gz

output_1=TU_sorted
output_1=AB_sorted

#Now run bwa mem alignment
bwa mem -t 32 ${REF_GENOME} ${TU_M_R1} ${TU_M_R2} -o TUBIGEN.sam

bwa mem -t 32 ${REF_GENOME} ${AB_F_R1} ${AB_F_R2} -o AB.sam

#Convert sam to bam (TU)
samtools view -b TUBIGEN.sam -t ${REF_GENOME}.fai -o TUBIGEN.bam

samtools sort --write-index -u TUBIGEN.bam -O BAM -o ${output_1}


#AB
samtools view -b AB.sam -t ${REF_GENOME}.fai -o AB.bam

samtools sort --write-index -u AB.bam -O BAM -o ${output_2}















