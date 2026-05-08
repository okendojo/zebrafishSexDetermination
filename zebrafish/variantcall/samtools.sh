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
module add bcftools

#move to the working directory

cd /data/okendojo/paradisfishProject/snv/fastq_files

# Define file paths
reference_genome="/data/okendojo/paradisfishProject/snv/assembly/macOpe2Assembly.fasta"  # Replace with the path to your reference genome
output_directory="samtoolsVariants"    # Directory to store output files


# Create the output directory if it doesn't exist
mkdir -p "$output_directory"

#samtools sort  macOpe_variant.bam -o aligned_sorted.bam
#samtools index aligned_sorted.bam


bcftools mpileup -Ov -f $reference_genome macOpe_variant_dedup.bam | bcftools call -mv -Ov -o macOpe2_variants.vcf

mv macOpe2_variants.vcf $output_directory
