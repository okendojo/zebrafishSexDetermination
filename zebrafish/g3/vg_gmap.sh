#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=100g
#SBATCH --ntasks-per-core=1
#SBATCH --time=90:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=gmap

#load modules

module load pangenome


cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/vcf_data

mkdir -p vgmapTempDir 

#Set file paths
rna_fastq1="/data/okendojo/zebrafish/data/g3/rna_sequences/trimmed/forwardreads.fq.gz"
rna_fasta2="/data/okendojo/zebrafish/data/g3/rna_sequences/trimmed/reversereads.fq.gz"

# Map simulated RNA-seq reads using vg mpmap
vg mpmap -n rna -t 4 -x rpvg_mpmap_run.spliced.xg -g rpvg_mpmap_run.spliced.gcsa -F GAMP -d rpvg_mpmap_run.spliced.dist -f $rna_fastq1 -f $rna_fasta2 > rpvg_mpmap_run.spliced.gamp
