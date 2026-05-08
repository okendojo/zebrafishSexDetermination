#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=300g
#SBATCH --ntasks-per-core=1
#SBATCH --time=90:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=vg_rpvg

#load modules

module load pangenome


cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/vcf_data

mkdir -p rpvgResults 

#Set file paths
ref="/data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.dna.primary_assembly.fa"
gtf="/data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.111.gtf"
vcf="/data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/vcf_data/timePointPhased.vcf.gz"



vg autoindex --threads 24 --workflow rpvg --verbosity 2 --prefix G3_rpvg_run --tmp-dir rpvgResults --ref-fasta ${ref} --vcf ${vcf} --tx-gff $gtf 


