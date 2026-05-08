#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=100g
#SBATCH --ntasks-per-core=1
#SBATCH --time=180:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=second_gmap2

#load modules

module load pangenome


cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/vcf_data


#Set file paths
fastq1="/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/bamFiles/merged_fastq.gz"


# Map simulated RNA-seq reads using vg mpmap
vg mpmap -n RNA -t 30 -x rpvg_mpmap_run.spliced.xg -g rpvg_mpmap_run.spliced.gcsa -F GAMP -d rpvg_mpmap_run.spliced.dist -i -f $fastq1 > rpvg_mpmap_run_3.spliced.gamp
