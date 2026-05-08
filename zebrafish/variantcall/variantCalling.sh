#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=64g
#SBATCH --ntasks-per-core=1
#SBATCH --time=96:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=variantCalling

#load the module
module load bwa
module load samtools
module load GATK
module load picard

#move to the working directory

cd /data/okendojo/zebrafish/data/fish11/variantCalling

# Define file paths
reference_genome="/data/Zebrafish_T2T/fish11/resolved_assembly/fish11.fasta"  # Replace with the path to your reference genome
output_directory="variant_output"    # Directory to store output files
read_1=/data/Zebrafish_T2T/fish11/illum/HM33TDSXY_19100322_S232_L004_R1_001.fastq.gz
read_2=/data/Zebrafish_T2T/fish11/illum/HM33TDSXY_19100322_S232_L004_R2_001.fastq.gz
sample2_1=/data/Zebrafish_T2T/fish11/illum/HTJLYDSXY_19100322_S16_L001_R1_001.fastq.gz
sample2_2=/data/Zebrafish_T2T/fish11/illum/HTJLYDSXY_19100322_S16_L001_R2_001.fastq.gz
sample3_1=/data/Zebrafish_T2T/fish11/illum/HTJLYDSXY_19100322_S16_L002_R1_001.fastq.gz
sample3_2=/data/Zebrafish_T2T/fish11/illum/HTJLYDSXY_19100322_S16_L002_R2_001.fastq.gz

# Create the output directory if it doesn't exist
mkdir -p "$output_directory"

# Step 1: Read Alignment with BWA
bwa mem -t 24 "$reference_genome" "$read_1" "$read_2" > "$output_directory/aligned_1.sam"

bwa mem -t 24 "$reference_genome" "$sample2_1" "$sample2_2" > "$output_directory/aligned_2.sam"

bwa mem -t 24 "$reference_genome" "$sample3_1" "$sample3_1" > "$output_directory/aligned_3.sam"

# Step 2: Convert SAM to BAM
samtools view -S -b "$output_directory/aligned_1.sam" > "$output_directory/aligned_1.bam"
samtools view -S -b "$output_directory/aligned_2.sam" > "$output_directory/aligned_2.bam"
samtools view -S -b "$output_directory/aligned_3.sam" > "$output_directory/aligned_3.bam"


# Step 3: Sort and index BAM file
samtools sort -o "$output_directory/aligned_1_sorted.bam" "$output_directory/aligned_1.bam"
samtools sort -o  "$output_directory/aligned_2_sorted.bam" "$output_directory/aligned_2.bam"
samtools sort -o  "$output_directory/aligned_3_sorted.bam" "$output_directory/aligned_3.bam" 

# step 3.1 Merge bam files
samtools merge -O BAM -@24 -o $output_directory/aligned_sorted.bam "$output_directory/*_sorted.bam"

#Add readgroups and index bam file

java -Xmx4g -XX:ParallelGCThreads=10 -jar $PICARDJARPATH/picard.jar AddOrReplaceReadGroups -I $output_directory/aligned_sorted.bam  -O $output_directory/aligned_sorted_RG.bam -RGID 4 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM 20 -SO coordinate --CREATE_INDEX true

# Step 4: Mark Duplicates with Picard
java -Xmx4g -XX:ParallelGCThreads=10 -jar $PICARDJARPATH/picard.jar MarkDuplicates -I $output_directory/aligned_sorted_RG.bam -O $output_directory/aligned_marked.bam -M $output_directory/marked_dup_metrics.txt

# Step 5: Base Quality Score Recalibration (BQSR) with GATK
gatk BaseRecalibrator -R "$reference_genome" -I "$output_directory/aligned_marked.bam" --known-sites /data/okendojo/zebrafish/data/fish11/variantCalling/GRCz11_known_variants.vcf.gz -O "$output_directory/recal_data.table"

gatk ApplyBQSR -R "$reference_genome" -I "$output_directory/aligned_marked.bam" --bqsr "$output_directory/recal_data.table" -O "$output_directory/aligned_recalibrated.bam"

# Step 6: Variant Calling with GATK
gatk HaplotypeCaller -R "$reference_genome" -I "$output_directory/aligned_recalibrated.bam" -O "$output_directory/raw_variants.vcf"

# Step 7: Variant Filtering with GATK
gatk VariantFiltration -R "$reference_genome" -V "$output_directory/raw_variants.vcf" -O "$output_directory/filtered_variants.vcf" \
  --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
  --filter-name "Filter"

# Step 8: Generate a final VCF file
gatk SelectVariants -R "$reference_genome" -V "$output_directory/filtered_variants.vcf" -O "$output_directory/final_variants.vcf" --exclude-filtered

echo "Variant calling and filtering complete."

