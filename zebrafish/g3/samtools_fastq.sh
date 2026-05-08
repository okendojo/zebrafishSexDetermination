#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=200g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=bamTofastq

#load modules
module add samtools

cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/bamFiles

for bam in *.bam ; do samtools fastq -1 ${bam}_R1_fastq.gz -2 ${bam}_R2_fastq.gz -n -@ 24 ${bam} ; done 
