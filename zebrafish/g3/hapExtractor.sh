#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=bcftoolsExtract



#Load the module
module load deeptools
module load samtools
module add bcftools

cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants 

REF="/data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.dna.primary_assembly.fa"

bcftools consensus -f ${REF} AB.vcf.gz -o AB.fasta -H 1


bcftools consensus -f ${REF} TL.vcf.gz -o TL.fasta -H 1


bcftools consensus -f ${REF} TU.vcf.gz -o TU.fasta -H 1


bcftools consensus -f ${REF} WIK.vcf.gz -o WIK.fasta -H 1


