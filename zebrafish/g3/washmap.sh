#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=90:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=what2

#load modules
module load samtools
module load whatshap

cd /data/okendojo/zebrafish/data/g3/eQTL/sarek_variants/sarekvariants/variant_calling/deepvariant/combinedVCFs

ref="/data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.dna.primary_assembly.fa"
AB_vcf="/data/okendojo/zebrafish/data/g3/eQTL/sarek_variants/sarekvariants/variant_calling/deepvariant/AB_strain/AB_strain.deepvariant.vcf.gz"
TU_vcf="/data/okendojo/zebrafish/data/g3/eQTL/sarek_variants/sarekvariants/variant_calling/deepvariant/TU_strain/TU_strain.deepvariant.vcf.gz"
TL_vcf="/data/okendojo/zebrafish/data/g3/eQTL/sarek_variants/sarekvariants/variant_calling/deepvariant/TL_strain/TL_strain.deepvariant.vcf.gz"
WIK_vcf="/data/okendojo/zebrafish/data/g3/eQTL/sarek_variants/sarekvariants/variant_calling/deepvariant/WIK_strain/WIK_strain.deepvariant.vcf.gz"

AB_bam="/data/okendojo/zebrafish/data/g3/eQTL/sarek_variants/sarekvariants/preprocessing/markduplicates/AB_strain/AB_strain.md.cram"
TU_bam="/data/okendojo/zebrafish/data/g3/eQTL/sarek_variants/sarekvariants/preprocessing/markduplicates/TU_strain/TU_strain.md.cram"
TL_bam="/data/okendojo/zebrafish/data/g3/eQTL/sarek_variants/sarekvariants/preprocessing/markduplicates/TL_strain/TL_strain.md.cram"
WIK_bam="/data/okendojo/zebrafish/data/g3/eQTL/sarek_variants/sarekvariants/preprocessing/markduplicates/WIK_strain/WIK_strain.md.cram"


whatshap phase --output f1_phasedsec.vcf.gz --output-read-list readsUsedphase.txt --reference ${ref} ${AB_vcf} $TU_vcf $TL_vcf $WIK_vcf $AB_bam $TU_bam $TL_bam $WIK_bam

