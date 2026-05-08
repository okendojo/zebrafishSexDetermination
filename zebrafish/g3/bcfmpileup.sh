#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=120g
#SBATCH --ntasks-per-core=1
#SBATCH --time=96:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=RNA_variantCall

#load module
module load bcftools

cd /data/okendojo/zebrafish/data/g3/illumina/WIK_M_CB/qc

Ref_FASTA="/data/okendojo/zebrafish/data/g3/illumina/WIK_M_CB/qc/GRCz11_genomic.fasta"
#Ref_FASTA="/data/Zebrafish_T2T/fish11/resolved_assembly/NHGRI_danRer_Tu_11.1.fasta"
BAM="/data/okendojo/zebrafish/data/g3/illumina/WIK_M_CB/qc/bamList.txt"

# 1, call variants
bcftools mpileup -Ou --redo-BAQ --min-BQ 20  --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f ${Ref_FASTA} -b ${BAM} | bcftools call -mv -Ov -o WIK.vcf

#SNP-Based Genome Binning
#MetaBAT will output bins (sets of contigs/scaffolds) that represent individual genomes based on the SNP information.

#metabat -i path/to/your/genome_assembly.fasta -a rnaVariants.vcf -o bins_dir

#Filter and normalize variants
#gatk SelectVariants -R reference_genome.fasta -V raw_variants.vcf -select-type SNP -O snps.vcf
#bcftools norm -m - -w 3 -o snps_normalized.vcf snps.vcf


#Genome binning
#This command creates windows of 10,000 base pairs with a step size of 5,000, and then counts the number of SNPs in each window.

#bedtools makewindows -g reference_genome.fasta.fai -w 10000 -s 5000 | bedtools map -a - -b snps_normalized.vcf -c 4 -o count > snp_density.bed


#visualize SNPS
# Example R code to plot SNP density
#snp_data <- read.table("snp_density.bed", header = FALSE, col.names = c("chrom", "start", "end", "snp_count"))
#plot(snp_data$chrom, snp_data$snp_count, type = "l", xlab = "Genomic Position", ylab = "SNP Count")


