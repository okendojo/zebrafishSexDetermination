#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=90:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=vg_index

#load modules

module load pangenome

#define file paths
cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/vcf_data

mkdir -p vgIndex

vg index --temp-dir vgIndex -t 24 -p -v timePointPhased.vcf.gz -g timePointPhased.gcsa -G f2.gbwt 
