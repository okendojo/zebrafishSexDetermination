#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=300g
#SBATCH --ntasks-per-core=1
#SBATCH --time=90:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=tmp_map

#load modules

module load pangenome


cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/vcf_data

mkdir -p tmp_vg_giraffeTempDir 

#Set file paths
ref="/data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.dna.primary_assembly.fa"
gtf="/data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.111.gtf"
vcf="/data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/vcf_data/timePointPhased.vcf.gz"



vg autoindex --target-mem 250G --threads 24 --workflow map --verbosity 2 --prefix tmp_G3_giraffe_run --tmp-dir tmp_vg_giraffeTempDir --ref-fasta ${ref} --vcf ${vcf} --tx-gff $gtf 


