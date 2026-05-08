

#Align DNA reads to the reference genome:
bwa mem -t 8 reference_genome.fa DNA_reads.fq | samtools sort -o DNA_reads.bam -

#Align RNA-seq reads to the reference genome:
STAR --genomeDir /path/to/genome/index --readFilesIn RNAseq_reads.fq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix RNAseq_reads

#Convert DNA BAM to BED format:
bedtools bamtobed -i DNA_reads.bam > DNA_reads.bed

#Intersect DNA and RNA-seq reads:
bedtools intersect -a DNA_reads.bed -b RNAseq_reads.bam -wa -wb > intersected_reads.bed

#The key idea is to use the mapping positions of DNA and RNA-seq reads to associate DNA reads with RNA transcripts or genomic regions


#!/bin/bash

# Paths to the VCF files (replace these with your actual file paths)
vcf1="file1.vcf"
vcf2="file2.vcf"
vcf3="file3.vcf"
vcf4="file4.vcf"

# Output directory
output_dir="variant_comparison_output"
mkdir -p "$output_dir"

# Compare VCF files and extract unique variants
bcftools isec -n=4 -c all \
  -o "$output_dir/unique_variants_all.vcf" \
  "$vcf1" "$vcf2" "$vcf3" "$vcf4"

# Extract unique variants from each file
for i in {1..4}; do
  bcftools view -Ov -o "$output_dir/unique_variants_file${i}.vcf" \
    "$output_dir/0000.vcf" "$output_dir/000${i}.vcf"
done

echo "Comparison completed. Unique variants for each file are in $output_dir."

