#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=300g
#SBATCH --ntasks-per-core=1
#SBATCH --time=90:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=gbwt

#load modules

module load pangenome


cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/vcf_data

mkdir -p vgmapTempDir 

#Set file paths
ref="/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/vcf_data/chr1_25.fasta"
gtf="/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/vcf_data/genes.gtf"
vcf="/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/vcf_data/chr1_25.vcf.gz"

vg gbwt -x  g3_mpmap_run.spliced.vg -v $vcf -o g3_mpmap_run.spliced.gbwt 
