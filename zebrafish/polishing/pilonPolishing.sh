#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=pilonPolishing

#Load modules
module add samtools
module add bwa
module add pilon

cd /data/okendojo/zebrafish/data/fish11/polishing

# Define file paths
genome_assembly="/data/Zebrafish_T2T/fish11/resolved_assembly/fish11.fasta"  # Replace with the path to your original genome assembly
short_reads="/data/Zebrafish_T2T/fish11/illum/illumina.fastq.gz"           # Replace with the path to your short-read data (e.g., Illumina reads)
output_directory="polished_assembly"       # Directory to store polished assembly

# Create the output directory if it doesn't exist
##mkdir -p "$output_directory"

# Step 1: Index the original assembly (if not already indexed)
##bwa index "$genome_assembly"

# Step 2: Align short reads to the original assembly
##bwa mem -t 24 "$genome_assembly" "$short_reads" > "$output_directory/aligned_reads.sam"

# Step 3: Convert SAM to BAM
##samtools view -S -b "$output_directory/aligned_reads.sam" > "$output_directory/aligned_reads.bam"

# Step 4: Sort and index the BAM file
##samtools sort -o "$output_directory/aligned_reads_sorted.bam" "$output_directory/aligned_reads.bam"
##samtools index "$output_directory/aligned_reads_sorted.bam"

# Step 5: Use Pilon for assembly polishing
pilon --genome /data/Zebrafish_T2T/fish11/resolved_assembly/fish11.fasta --frags /data/okendojo/zebrafish/data/fish11/polishing/polished_assembly/aligned_reads_sorted.bam --output /data/okendojo/zebrafish/data/fish11/polishing/polished_assembly --threads 24 --vcf --tracks 

echo "Genome assembly polishing complete."

