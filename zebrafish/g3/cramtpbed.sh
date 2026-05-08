#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=200g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=cramTo:bam


#module loadings
module load samtools
module load bedtools

cd /data/okendojo/zebrafish/data/g3/eQTL/sarek_variants/dna_variants/variant_calling/deepvariant 

samtools view -b -T  /data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.dna.primary_assembly.fa /data/okendojo/zebrafish/data/g3/eQTL/sarek_variants/dna_variants/preprocessing/markduplicates/AB_strain/AB_strain.md.cram -o AB.bam

samtools view -b -T  /data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.dna.primary_assembly.fa  /data/okendojo/zebrafish/data/g3/eQTL/sarek_variants/dna_variants/preprocessing/markduplicates/TL_strain/TL_strain.md.cram -o TL.bam

samtools view -b -T  /data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.dna.primary_assembly.fa  /data/okendojo/zebrafish/data/g3/eQTL/sarek_variants/dna_variants/preprocessing/markduplicates/TU_strain/TU_strain.md.cram -o TU.bam
 
samtools view -b -T  /data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.dna.primary_assembly.fa  /data/okendojo/zebrafish/data/g3/eQTL/sarek_variants/dna_variants/preprocessing/markduplicates/WIK_strain/WIK_strain.md.cram -o WIK.bam

zcat AB_strain/AB_strain.deepvariant.vcf.gz | awk '{if($1 !~ /^#/) print $1"\t"$2-1"\t"$2"\t"$3}' > AB.bed
zcat TL_strain/TL_strain.deepvariant.vcf.gz | awk '{if($1 !~ /^#/) print $1"\t"$2-1"\t"$2"\t"$3}' > TL.bed
zcat TU_strain/TU_strain.deepvariant.vcf.gz  |  awk '{if($1 !~ /^#/) print $1"\t"$2-1"\t"$2"\t"$3}' > TU.bed
zcat WIK_strain/WIK_strain.deepvariant.vcf.gz | awk '{if($1 !~ /^#/) print $1"\t"$2-1"\t"$2"\t"$3}' > WIK.bed


bedtools intersect -c -a AB.bed -b AB.bam > AB_read_counts.txt
bedtools intersect -c -a TU.bed -b TU.bam > TU_read_counts.txt
bedtools intersect -c -a TL.bed -b TL.bam > TL_read_counts.txt
bedtools intersect -c -a WIK.bed -b WIK.bam > WIK_read_counts.txt


