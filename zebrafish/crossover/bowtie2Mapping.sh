#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=bowtie2_ZF

cd /data/okendojo/zebrafish/data/G3/TU_AB_alignment/treamedReads

#Load the modules
module load bowtie
module load samtools

#See the data paths
REF_GENOME=/data/okendojo/zebrafish/refGenome/GRCz11_genomic.fasta
INDEX=/data/okendojo/zebrafish/refGenome/index/index
TU_M_R1=/data/okendojo/zebrafish/data/G3/TU_M_CB/qc_results/HFVYJDSX5_19412541_S50_L002_R1_001_trimmed.fq.gz
TU_M_R2=/data/okendojo/zebrafish/data/G3/TU_M_CB/qc_results/HFVYJDSX5_19412541_S50_L002_R2_001_trimmed.fq.gz
AB_F_R1=/data/okendojo/zebrafish/data/G3/AB_F_CB/qc_results/HFVYJDSX5_19412535_S54_L002_R1_001_trimmed.fq.gz
AB_F_R2=/data/okendojo/zebrafish/data/G3/AB_F_CB/qc_results/HFVYJDSX5_19412535_S54_L002_R2_001_trimmed.fq.gz

#Run bowtie2

bowtie2 --no-unal -p 32 -x ${INDEX} -1 ${TU_M_R1} -2 ${TU_M_R2} -S TU.sam

bowtie2 --no-unal -p 32 -x ${INDEX} -1 ${AB_F_R1} -2 ${AB_F_R2} -S AB.sam


#Convert sam to bam
samtools view -b -u -o TU.bam -t ${REF_GENOME}.fai TU.sam
samtools view -b -u -o AB.bam -t ${REF_GENOME}.fai AB.sam

#Run stampy to improve the alignment of the indels
#Set the path variables for the genome index and the harsh table
index=/data/okendojo/zebrafish/refGenome/index/genome
hashtable=/data/okendojo/zebrafish/refGenome/index/genome


#Run sampy to select good reads only

stampy.py -g ${index} -h ${hashtable} -t 32 --bamkeepgoodreads -f sam -o TU_stampyCleaned.sam -M TU.bam

stampy.py -g ${index} -h ${hashtable} -t 32 --bamkeepgoodreads -f sam -o AB_stampCleaned.sam -M AB.bam
 
#remove the duplicate reads from the bam files

echo "........we are now starting to remove the duplicate reads at " `date`
for file in *_stampCleaned.sam ;

do

samtools rmdup -sS  --output-fmt BAM --verbosity 4 ${file} ${file}_dedup.bam

done

echo "........duplicate removal is done at" `date`


#Sort bam files

samtools sort -u -m 2G -o TU_final.bam  --threads 20 -O BAM TU_stampyCleaned.sam
samtools sort -u -m 2G -o AB_final.bam  --threads 20 -O BAM AB_stampCleaned.sam










