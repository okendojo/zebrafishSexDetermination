#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=230g
#SBATCH --ntasks-per-core=1
#SBATCH --time=140:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=mapIsoSeq


#Load the modules
ml minimap2 samtools

cd /data/okendojo/zebrafish/data/g3/sex_project/asm/t2tMap/isoseq

minimap2 -ax splice:hq -uf --secondary=no -t 24 /data/okendojo/zebrafish/data/g3/sex_project/asm/t2tMap/NHGRI_Fish6.fasta /vf/users/okendojo/zebrafish/data/g3/iso-seq/sra_isoseq/m84270_240911_150713_s1.skera.flnc.fastq.gz | samtools sort -@ 24 -o testise.bam

samtools index testise.bam


minimap2 -ax splice:hq -uf --secondary=no -t 24 /data/okendojo/zebrafish/data/g3/sex_project/asm/t2tMap/NHGRI_Fish6.fasta /vf/users/okendojo/zebrafish/data/g3/iso-seq/sra_isoseq/m84270_240911_210533_s4.skera.flnc.fastq.gz | samtools sort -@ 24 -o overy.bam

samtools index ovary.bam
