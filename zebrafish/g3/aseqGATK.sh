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

cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/ASEQ


fasta="/data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.dna.primary_assembly.fa"
AB_T0=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/preprocessing/AB_T0/AB_T0.recal.bam
AB_T12=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/preprocessing/AB_T12/AB_T12.recal.bam
AB_T24=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/preprocessing/AB_T24/AB_T24.recal.bam
AB_T6=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/preprocessing/AB_T6/AB_T6.recal.bam

TL_T0=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/preprocessing/TL_T0/TL_T0.recal.bam
TL_T12=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/preprocessing/TL_T12/TL_T12.recal.bam
TL_T24=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/preprocessing/TL_T24/TL_T24.recal.bam
TL_T6=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/preprocessing/TL_T6/TL_T6.recal.bam

TU_T0=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/preprocessing/TU_T0/TU_T0.recal.bam
TU_T12=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/preprocessing/TU_T12/TU_T12.recal.bam
TU_T24=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/preprocessing/TU_T24/TU_T24.recal.bam
TU_T6=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/preprocessing/TU_T6/TU_T6.recal.bam

WIK_T0=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/preprocessing/WIK_T0/WIK_T0.recal.bam
WIK_T12=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/preprocessing/WIK_T12/WIK_T12.recal.bam
WIK_T24=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/preprocessing/WIK_T24/WIK_T24.recal.bam
WIK_T6=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/preprocessing/WIK_T6/WIK_T6.recal.bam

#VCFs
AB_T0v=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/variant_calling/AB_T0/AB_SNPs_out.vcf.gz
AB_T12v=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/variant_calling/AB_T12/out.vcf.gz
AB_T24v=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/variant_calling/AB_T24/out.vcf.gz
AB_T6v=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/variant_calling/AB_T6/out.vcf.gz


TL_T0v=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/variant_calling/TL_T0/out.vcf.gz
TL_T12v=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/variant_calling/TL_T12/out.vcf.gz
TL_T24v=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/variant_calling/TL_T24/out.vcf.gz
TL_T6v=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/variant_calling/TL_T6/out.vcf.gz

TU_T0v=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/variant_calling/TU_T0/out.vcf.gz
TU_T12v=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/variant_calling/TU_T12/out.vcf.gz
TU_T24v=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/variant_calling/TU_T24/out.vcf.gz
TU_T6v=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/variant_calling/TU_T6/out.vcf.gz

WIK_T0v=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/variant_calling/WIK_T0/out.vcf.gz
WIK_T12v=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/variant_calling/WIK_T12/out.vcf.gz
WIK_T24v=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/variant_calling/WIK_T24/out.vcf.gz
WIK_T6v=/vf/users/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/rnavarAnalysis/rna_variants/variant_calling/WIK_T6/out.vcf.gz




gatk ASEReadCounter -R $fasta -I ${AB_T0} --variant ${AB_T0v}  --min-depth-of-non-filtered-base 10 --min-base-quality 10 --min-mapping-quality 10  --output-format CSV --output AB_T0_ASEresults.csv 

gatk ASEReadCounter -R $fasta -I ${AB_T12} --variant ${AB_T12v} --min-base-quality 10 --min-depth-of-non-filtered-base 10 --min-mapping-quality 10    --output-format CSV --output AB_T12_ASEresults.csv

gatk ASEReadCounter -R $fasta -I ${AB_T24} --variant ${AB_T24v} --min-base-quality 10 --min-mapping-quality 10 --min-depth-of-non-filtered-base 10  --output-format CSV --output AB_T24_ASEresults.csv

gatk ASEReadCounter -R $fasta -I ${AB_T6} --variant ${AB_T6v} --min-base-quality 10 --min-mapping-quality 10 --min-depth-of-non-filtered-base 10   --output-format CSV --output AB_T6_ASEresults.csv

gatk ASEReadCounter -R $fasta -I ${TL_T0} --variant ${TL_T0v} --min-base-quality 10 --min-mapping-quality 10 --min-depth-of-non-filtered-base 10  --output-format CSV --output TL_T0_ASEresults.csv

gatk ASEReadCounter -R $fasta -I ${TL_T12} --variant ${TL_T12v} --min-base-quality 10 --min-mapping-quality 10 --min-depth-of-non-filtered-base 10   --output-format CSV --output TL_T12_ASEresults.csv

gatk ASEReadCounter -R $fasta -I ${TL_T24} --variant ${TL_T24v} --min-base-quality 10 --min-mapping-quality 10 --min-depth-of-non-filtered-base 10   --output-format CSV --output TL_T24_ASEresults.csv

gatk ASEReadCounter -R $fasta -I ${TL_T6} --variant ${TL_T6v} --min-base-quality 10 --min-mapping-quality 10 --min-depth-of-non-filtered-base 10  --output-format CSV  --output TL_T6_ASEresults.csv

gatk ASEReadCounter -R $fasta -I ${TU_T0} --variant ${TU_T0v} --min-base-quality 10 --min-mapping-quality 10 --min-depth-of-non-filtered-base 10   --output-format CSV --output TU_T0_ASEresults.csv

gatk ASEReadCounter -R $fasta -I ${TU_T12} --variant ${TU_T12v} --min-base-quality 10 --min-mapping-quality 10 --min-depth-of-non-filtered-base 10  --output-format CSV --output TU_T12_ASEresults.csv

gatk ASEReadCounter -R $fasta -I ${TU_T24} --variant ${TU_T24v} --min-base-quality 10 --min-mapping-quality 10 --min-depth-of-non-filtered-base 10   --output-format CSV --output TU_T24_ASEresults.csv

gatk ASEReadCounter -R $fasta -I ${TU_T6} --variant ${TU_T6v} --min-base-quality 10 --min-mapping-quality 10 --min-depth-of-non-filtered-base 10  --output-format CSV --output TU_T6_ASEresults.csv

gatk ASEReadCounter -R $fasta -I ${WIK_T0} --variant ${WIK_T0v} --min-base-quality 10 --min-mapping-quality 10 --min-depth-of-non-filtered-base 10  --output-format CSV --output WIK_T0_ASEresults.csv

gatk ASEReadCounter -R $fasta -I ${WIK_T12} --variant ${WIK_T12v} --min-base-quality 10 --min-mapping-quality 10 --min-depth-of-non-filtered-base 10  --output-format CSV --output WIK_T12_ASEresults.csv

gatk ASEReadCounter -R $fasta -I ${WIK_T24} --variant ${WIK_T24v} --min-base-quality 10 --min-mapping-quality 10 --min-depth-of-non-filtered-base 10 --output-format CSV   --output WIK_T24_ASEresults.csv

gatk ASEReadCounter -R $fasta -I ${WIK_T6} --variant ${WIK_T6v} --min-base-quality 10 --min-mapping-quality 10 --min-depth-of-non-filtered-base 10   --output-format CSV --output WIK_T6_ASEresults.csv


