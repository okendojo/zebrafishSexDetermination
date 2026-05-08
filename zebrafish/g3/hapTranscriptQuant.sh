#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=90:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=haptransQuant

#load modules

module load pangenome


cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/vcf_data

script="/home/okendojo/scripts/zebrafish/g3/PDTG-project/scripts/haplotype-transcript"
ref="/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/vcf_data/chr1_25.fasta"
gtf="/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/vcf_data/genes.gtf"
phasedVCF="/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/vcf_data/chr1_25.vcf.gz"
fq1="/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/bamFiles/TL_T12_variants.bam_R1_fastq.gz"
fq2="/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/bamFiles/TL_T12_variants.bam_R2_fastq.gz"

#rna_fastq1="/data/okendojo/zebrafish/data/g3/rna_sequences/trimmed/*_R1_001_trimmed.fq.gz"
#rna_fasta2="/data/okendojo/zebrafish/data/g3/rna_sequences/trimmed/*_R2_001_trimmed.fq.gz"
G3_rna_graph="/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/vcf_data/rpvg_mpmap_run.spliced.xg"


sh ${script}/Mainscript.autoindexRun.sh chr1.fasta chr1.vcf.gz genes_chr1.gtf rpvg_mpmap_chr1.spliced rpvg_mpmap_chr1.spliced run2Hap/ /vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/bamFiles/TL_T12_variants.bam_R1_fastq.gz /vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/bamFiles/TL_T12_variants.bam_R2_fastq.gz
