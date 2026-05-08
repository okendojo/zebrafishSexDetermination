#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=120g
#SBATCH --ntasks-per-core=1
#SBATCH --time=96:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov

start=`date +%s`
echo $HOSTNAME

#load the module
module load bwa
module load samtools
module load dragmap/1.2.1

cd /data/okendojo/zebrafish/data/g3/illumina/bamFiles

REF=/data/okendojo/zebrafish/refGenome/GRCz11.fasta
AB=/data/okendojo/zebrafish/data/g3/illumina/AB_F_CB/raw
TL=/data/okendojo/zebrafish/data/g3/illumina/TL_F_CB/raw
TU=/data/okendojo/zebrafish/data/g3/illumina/TU_M_CB/raw
WIK=/data/okendojo/zebrafish/data/g3/illumina/WIK_M_CB/raw

bwa mem -t 32 ${REF} $AB/HFVYJDSX5_19412535_S54_L002_R1_001_trimmed.fq.gz $AB/HFVYJDSX5_19412535_S54_L002_R2_001_trimmed.fq.gz -o AB_F.sam

bwa mem -t 32 ${REF} $TL/HFVYJDSX5_19412527_S51_L002_R1_001_trimmed.fq.gz $TL/HFVYJDSX5_19412527_S51_L002_R2_001_trimmed.fq.gz -o TL_F.sam

bwa mem -t 32 $REF $TU/HFVYJDSX5_19412541_S50_L002_R1_001_trimmed.fq.gz $TU/HFVYJDSX5_19412541_S50_L002_R2_001_trimmed.fq.gz -o TU_M.sam

bwa mem -t 32 $REF $WIK/HFVYJDSX5_19412539_S55_L002_R1_001_trimmed.fq.gz $WIK/HFVYJDSX5_19412539_S55_L002_R2_001_trimmed.fq.gz -o WIK_M.sam

echo "samtools analysis is starting........."

for sam in *.sam
do

samtools view -S -b $sam -o ${sam}.bam

done

echo "Done converting sam to bam!"

for bam in *.bam

do

samtools sort $bam -o ${sam}_sorted.bam 

done

echo "Done sorting bam files"

for sortedbam *.sorted.bam

do

samtools index ${sortedbam}

done

echo "Done indexing the files"





