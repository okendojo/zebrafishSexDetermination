#!/bin/bash

# Define input files
reference_genome="reference.fasta"
bam_file="mapped_reads.bam"
filtered_variants_vcf="filtered_variants.vcf"

# Step 1: Create pileup file using samtools mpileup
samtools mpileup -f "$reference_genome" "$bam_file" > pileup.txt

# Step 2: Process pileup file to get allelic counts for filtered variants
awk 'BEGIN {OFS="\t"} {print $1,$2-1,$2}' "$filtered_variants_vcf" | \
    bedtools intersect -a - -b pileup.txt -wa -wb | \
    awk '{print $4,$5,$6,$7,$8,$9,$10}' > allelic_counts.txt

# Output format: Chromosome, Position, Reference base, Allelic depth for reference allele, Allelic depth for alternate allele(s)

# Optionally, filter allelic counts based on quality thresholds or other criteria
# Further analysis and interpretation of allelic counts

