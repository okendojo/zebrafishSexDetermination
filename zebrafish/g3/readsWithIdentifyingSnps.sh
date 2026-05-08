#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=90:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=readsWithSnps

#This script is used to identify how many RNA-seq reads align to transcripts that have identifying SNPs
#move to the right directory
cd /data/okendojo/zebrafish/data/g3/rna_sequences/gene_quantification/SNP_analysis

#Load the modules
module load GATK
module load samtools
module load subread 


# Define paths to input files
reference_genome="/data/okendojo/zebrafish/data/g3/rna_sequences/gene_quantification/ncbi_dataset/Danio_rerio.GRCz11.dna.primary_assembly.fa"
transcript_annotations="/data/okendojo/zebrafish/data/g3/rna_sequences/gene_quantification/ncbi_dataset/Danio_rerio.GRCz11.110.chr.gtf"
input_bam="/data/okendojo/zebrafish/data/g3/rna_sequences/gene_quantification/rnaseq_result/star_salmon/time_0_1.markdup.sorted.bam"  # BAM file with aligned RNA-seq reads

# Output directory
output_dir="output_dir"
mkdir -p $output_dir

# Step 1: Index the reference genome
#samtools faidx $reference_genome

# Step 2: Create a sequence dictionary for the reference genome
#gatk CreateSequenceDictionary -R $reference_genome -O $output_dir/reference.dict

# Step 3: Identify SNPs using GATK HaplotypeCaller
gatk HaplotypeCaller -R $reference_genome -I $input_bam -O $output_dir/raw_variants.vcf

# Step 4: Filter SNPs based on quality and coverage
gatk VariantFiltration -R $reference_genome -V $output_dir/raw_variants.vcf -O $output_dir/filtered_variants.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" --filter-name "snp_filter"

# Step 5: Extract positions of filtered SNPs
awk '/^#/ {print; next} $7 == "snp_filter" {print $1, $2}' $output_dir/filtered_variants.vcf > $output_dir/snp_positions.bed

# Step 6: Extract reads overlapping SNP positions
samtools view -L $output_dir/snp_positions.bed $input_bam -o $output_dir/snp_overlapping_reads.bam -O BAM ---reference $reference_genome 

# Step 7: Count reads per transcript
featureCounts -a $transcript_annotations -o $output_dir/read_counts_per_transcript.txt -F GTF $output_dir/snp_overlapping_reads.bam
