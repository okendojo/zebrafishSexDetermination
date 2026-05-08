#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=fish11_2_t2t


#Load the modules
module load bcftools
module load sniffles
module load whatshap
module load minimap2 samtools


#move to the working directory
cd /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/haplotypevariation

# Define input files
REFERENCE="/data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/asm_polished/NHGRI_Fish6_cons.fasta"   # Reference genome
READS="/data/okendojo/zebrafish/data/fish6/hifi/m54313U_220817_024630.hifi_reads.fastq.gz"             # PacBio/ONT long reads
BAM="aligned_reads.bam"              # Output BAM file
VCF="variants.vcf"                   # Variant Call File
SV_VCF="structural_variants.vcf"     # Structural Variants VCF
OUTPUT_PREFIX="zebrafish_analysis"

# Number of threads
THREADS=24

echo "Starting Zebrafish Haplotype & Structural Variation Analysis"

# Step 1: Index the Reference Genome
echo "Indexing the Reference Genome..."
samtools faidx $REFERENCE

# Step 2: Align Reads to Reference Genome using Minimap2
echo "Aligning reads using Minimap2..."
minimap2 -ax map-pb $REFERENCE $READS | samtools sort -o $BAM
samtools index $BAM

# Step 3: Call Variants (SNPs and Indels) using BCFtools
echo "Calling Variants using BCFtools..."
bcftools mpileup -Ou -f $REFERENCE $BAM | bcftools call -mv -Oz -o $VCF
bcftools index $VCF

# Step 4: Phase Variants using Whatshap
echo "Phasing Variants with Whatshap..."
whatshap phase --reference $REFERENCE -o phased_variants.vcf $VCF $BAM

# Step 5: Detect Structural Variations with Sniffles
echo "Detecting Structural Variations with Sniffles..."
sniffles --input $BAM --vcf $SV_VCF --reference $REFERENCE --threads $THREADS

# Step 6: Filter High-Quality Structural Variants
echo "Filtering high-confidence Structural Variants..."
bcftools filter -i 'QUAL>30' $SV_VCF > high_quality_SVs.vcf

# Step 7: Generate Summary Statistics
echo "Generating Summary Statistics..."
bcftools stats high_quality_SVs.vcf > SV_summary.txt

# Step 8: Prepare Data for Circos Visualization
echo "Formatting Structural Variants for Circos plot..."
grep -v "^#" high_quality_SVs.vcf | awk '{print $1"\t"$2"\t"$2+1"\t"$5}' > SV_circos.txt

echo "Haplotype & Structural Variation Analysis Completed!"

