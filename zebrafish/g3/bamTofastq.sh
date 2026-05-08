#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=120g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=bcftools_iseq

#load module
module load bcftools
module load subread
module load samtools

cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/bamFiles 

for bam in *.bam ; do samtools sort -n -o ${bam}_sorted.bam $bam ; done # sort bam files

for file in *_sorted.bam ; do samtools index $file ; done # index the sorted bam files

for index in *_sorted.bam ; do bedtools bamtofastq -i ${index}  -fq ${index}_R1.fq -fq2 ${index}_R2.fq ; done # convert bed to fastq

