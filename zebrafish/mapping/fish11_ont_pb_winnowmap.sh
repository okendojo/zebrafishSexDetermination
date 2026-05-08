#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=winnowmap

#load the modules

module load winnowmap
module load samtools
module load meryl

cd /data/okendojo/zebrafish/data/fish11/mapping

assembly=/vf/users/okendojo/zebrafish/data/fish11/gapfilling/fish11.fasta
ont=/vf/users/okendojo/zebrafish/data/fish11/ont/fastq_pass/concatenatedONT/ont.fasta
hifi=/vf/users/okendojo/zebrafish/data/fish11/hifi/hifi.fasta 

#===========================ONT mapping =============================
#Run the ONT mapping
winnowmap -t 24 -ax map-ont ${assembly} ${ont} > fish11_ontMapping.sam

echo "....now running samtools on ont assembly map"

samtools view -S -b fish11_ontMapping.sam -o fish11_ontMapping.bam
samtools sort fish11_ontMapping.bam -o fish11_ontMapping_sorted.bam
samtools index fish11_ontMapping_sorted.bam

#=============================PacBio=================================
#mapping PacBio hifi data

echo "....now mapping PacBio data"
winnowmap -t 24 -ax map-pb ${assembly} ${hifi} > fish11_pbMapping.sam


echo "Now running samtools.........."

samtools view -S -b fish11_pbMapping.sam -o fish11_pbMapping.bam
samtools sort fish11_pbMapping.bam -o fish11_pbMapping_sorted.bam
samtools index fish11_pbMapping_sorted.bam




















