#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=minimap2

#load module
module add minimap2
module load samtools
module load bowtie2 

cd /data/okendojo/zebrafish/data/g3/mapping

# 1. Build the assembly index
asm="/data/okendojo/zebrafish/data/g3/assembly/verkko_asm2/assembly.fasta"
bowtie2-build  -f ${asm} --threads 24

# 2. Reads alignment

#Set data paths
REF="/data/okendojo/zebrafish/data/g3/assembly/verkko_asm2/assembly.fasta"
T2T="/data/okendojo/zebrafish/data/g3/pacBio/*.fastq.gz"

minimap2 -x map-hifi -a -t 24  -N 1 --secondary=no $REF $T2T  | samtools view -@ 20 -bS - | samtools sort -@ 20 -o HiFi_ManualCurated.bam -

samtools view -h -F 3840 -bS HiFi_ManualCurated.bam | samtools sort -o HiFi_ManualCurated_clean.bam -
samtools depth -Q 1 -a HiFi_ManualCurated_clean.bam > HiFi_ManualCurated_depth_clean.txt
