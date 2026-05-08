# Assuming you have paired-end RNA-seq data in fastq format

# Indexing reference genome with HISAT2
hisat2-build reference_genome.fa reference_genome

# Aligning reads with HISAT2
hisat2 -x reference_genome -1 sample_R1.fastq -2 sample_R2.fastq -S sample.sam

# Converting SAM to BAM and sorting
samtools view -bS sample.sam | samtools sort -o sample_sorted.bam

# Using StringTie for transcript assembly
stringtie sample_sorted.bam -o transcripts.gtf -G reference_annotation.gtf

# Merge all samples if you have multiple time points
stringtie --merge -o merged_transcripts.gtf sample1.gtf sample2.gtf ...

# Example using DESeq2
library(DESeq2)

# Assuming you have a DESeqDataSet object named "dds"
dds <- DESeqDataSetFromHTSeqCount(sampleTable, directory = "path/to/htseq_counts", design = ~time_point)
dds <- DESeq(dds)
res <- results(dds)

# Assuming you have genomic data in BAM format
# Use GATK for variant calling
gatk HaplotypeCaller -R reference_genome.fa -I genomic_data.bam -O variants.vcf

# Assuming you have the DEGs in a file named "DEGs.txt"

# Annotate variants with genes
bedtools intersect -a variants.vcf -b reference_annotation.gtf -wa -wb > annotated_variants.bed

# Extract DEGs associated with variants
awk 'NR==FNR{a[$4]; next} $NF in a' <(cut -f9 annotated_variants.bed | grep -o 'gene_id "[^"]*' | sed 's/gene_id "//') DEGs.txt > DEGs_with_variants.txt

