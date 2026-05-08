#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=200g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=GATK_Annot

#load modules
module add GATK

cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/vcf_data

REF="/data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.dna.primary_assembly.fa"

gatk VariantAnnotator --variant filtred_T12.vcf.gz -R ${REF} --output  filtred_T12_annot.vcf.gz  -A SampleList 
