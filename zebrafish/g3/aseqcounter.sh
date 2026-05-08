#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=200g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=GATK_ASEcounter

#load modules
module add GATK

cd /data/okendojo/zebrafish/data/g3/eQTL/ASEReadCounter

fasta="/data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.dna.primary_assembly.fa"
BAM_T0="/data/okendojo/zebrafish/data/g3/eQTL/fullTimePoints/preprocessing/time_0/time_0.recal.bam"
BAM_T6="/data/okendojo/zebrafish/data/g3/eQTL/fullTimePoints/preprocessing/time_6/time_6.recal.bam"
BAM_T12="/data/okendojo/zebrafish/data/g3/eQTL/fullTimePoints/preprocessing/time_12/time_12.recal.bam"
BAM_T24="/data/okendojo/zebrafish/data/g3/eQTL/fullTimePoints/preprocessing/time_24/time_24.recal.bam"

VCF_T0="/data/okendojo/zebrafish/data/g3/eQTL/fullTimePoints/variant_calling/time_0/time_0_data.vcf.gz"
VCF_T6="/data/okendojo/zebrafish/data/g3/eQTL/fullTimePoints/variant_calling/time_6/time_6_variants.vcf.gz"
VCF_T12="/data/okendojo/zebrafish/data/g3/eQTL/fullTimePoints/variant_calling/time_12/time_12_variants.vcf.gz"
VCF_T24="/data/okendojo/zebrafish/data/g3/eQTL/fullTimePoints/variant_calling/time_24/time_24_variants.vcf.gz"

gatk ASEReadCounter -R $fasta -I ${BAM_T0} --variant ${VCF_T0}  --min-depth-of-non-filtered-base 10 --min-base-quality 10 --min-mapping-quality 20 --output-format CSV --output time_0.csv 

gatk ASEReadCounter -R $fasta -I ${BAM_T6} --variant ${VCF_T6} --min-base-quality 10 --min-depth-of-non-filtered-base 10 --min-mapping-quality 20  --output-format CSV  --output time_6.csv

gatk ASEReadCounter -R $fasta -I ${BAM_T12} --variant ${VCF_T12} --min-base-quality 10 --min-mapping-quality 20 --min-depth-of-non-filtered-base 10 --output-format CSV --output time_12.csv

gatk ASEReadCounter -R $fasta -I ${BAM_T24} --variant ${VCF_T24} --min-base-quality 10 --min-mapping-quality 20 --min-depth-of-non-filtered-base 10 --output-format CSV --output time_24.csv


