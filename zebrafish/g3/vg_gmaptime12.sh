#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=200g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=time12_gmap2

#load modules

module load pangenome


cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/vcf_data


#Set file paths
#fastq1="/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/bamFiles/merged_time12.fastq.gz"

fq1="/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/bamFiles/TL_T12_variants.bam_R1_fastq.gz"
fq2="/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/bamFiles/TL_T12_variants.bam_R2_fastq.gz"

# Map simulated RNA-seq reads using vg mpmap
vg mpmap -n RNA -t 30 -x rpvg_mpmap_chr1.spliced.xg -g rpvg_mpmap_chr1.spliced.gcsa -F GAMP -d rpvg_mpmap_chr1.spliced.dist -f $fq1 -f $fq2 > rpvg_mpmap_chr1.spliced.gamp
