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

cd /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/pur

minimap2 -a -x asm5 --eqx NHGRI_Fish6_cons.fasta GRCz11.fasta -o NHGRI_GRCz11.sorted.sam 
samtools view -Sb NHGRI_GRCz11.sorted.sam > NHGRI_GRCz11.sorted.bam
samtools sort -o NHGRI_GRCz11.sorted.bam
samtools index NHGRI_GRCz11.sorted.bam

