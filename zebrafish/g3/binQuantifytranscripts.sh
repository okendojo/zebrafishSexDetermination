# Define input files
bam_file="your_input.bam"
vcf_file="your_variants.vcf"
reference_genome="reference_genome.fa"

# Step 1: Sort and index the BAM file (if not already done)
samtools sort -o sorted.bam "${bam_file}"
samtools index sorted.bam

# Step 2: Create a BED file from the VCF file
# Assuming you are interested in coding regions only (modify as needed)
bcftools query -f '%CHROM\t%POS0\t%POS1\t%REF\t%ALT\n' "${vcf_file}" > variants.bed

# Step 3: Create a transcript BED file (modify as needed)
# Assuming you have a transcript BED file, or you can generate it based on your reference annotation
# Example: Use Gencode annotation to create a transcript BED file
# gencode_annotation.gtf -> Convert to BED format
# gtfToGenePred gencode_annotation.gtf gencode_annotation.genepred
# genePredToBed gencode_annotation.genepred gencode_transcripts.bed

# Step 4: Intersect the variant BED file with the transcript BED file
bedtools intersect -a variants.bed -b gencode_transcripts.bed -wa -wb > intersected_variants.bed

# Step 5: Extract reads overlapping with variant positions
bedtools intersect -a sorted.bam -b intersected_variants.bed > reads_overlapping_variants.bam

# Step 6: Quantify transcripts using GATK
gatk BaseRecalibrator -I reads_overlapping_variants.bam -R "${reference_genome}" -O recal_data.table --known-sites "${vcf_file}"
gatk ApplyBQSR -I reads_overlapping_variants.bam -R "${reference_genome}" --bqsr-recal-file recal_data.table -O recalibrated.bam

# Step 7: Quantify transcripts using featureCounts
# Assuming you have a transcript annotation file in GTF format (modify as needed)
featureCounts -a gencode_annotation.gtf -o counts.txt -R BAM recalibrated.bam

