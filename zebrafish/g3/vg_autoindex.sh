#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=90:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=T1_autoindex

#load modules

module load pangenome


cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/vcf_data

mkdir -p vgmapTempDir 

#Set file paths
ref="/data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.dna.primary_assembly.fa"
gtf="/data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.111.gtf"
vcf="/data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/vcf_data/timePointPhased.vcf.gz"



vg autoindex --threads 24 --workflow mpmap --verbosity 2 --prefix time_12vgp --tmp-dir vgmapTempDir --ref-fasta ${ref} --vcf ${vcf}  --tx-gff ${gtf}


