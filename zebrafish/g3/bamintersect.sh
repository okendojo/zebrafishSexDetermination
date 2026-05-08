#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=120g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=bcftools_iseq

#load module
module load bcftools
module load subread

cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis

TIME_0="/data/okendojo/zebrafish/data/g3/eQTL/fullTimePoints/preprocessing/time_0/time_0.recal.bam"
TIME_6="/data/okendojo/zebrafish/data/g3/eQTL/fullTimePoints/preprocessing/time_6/time_6.recal.bam"
TIME_12="/data/okendojo/zebrafish/data/g3/eQTL/fullTimePoints/preprocessing/time_12/time_12.recal.bam"
TIME_24="/data/okendojo/zebrafish/data/g3/eQTL/fullTimePoints/preprocessing/time_24/time_24.recal.bam"
GTF="/data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.111.gtf"


#Run time zero 
bedtools intersect -a ${TIME_0} -b AB_geneIntersect.bed > AB_T0_variants.bam

bedtools intersect -a ${TIME_0}  -b TU_geneIntersect.bed > TU_T0_variants.bam

bedtools intersect -a ${TIME_0} -b TL_geneIntersect.bed > TL_T0_variants.bam

bedtools intersect -a  ${TIME_0} -b WIK_geneIntersect.bed > WIK_T0_variants.bam



bedtools intersect -a  $TIME_6 -b AB_geneIntersect.bed > AB_T6_variants.bam

bedtools intersect -a  $TIME_6 -b TU_geneIntersect.bed > TU_T6_variants.bam

bedtools intersect -a $TIME_6 -b TL_geneIntersect.bed > TL_T6_variants.bam

bedtools intersect -a $TIME_6 -b WIK_geneIntersect.bed > WIK_T6_variants.bam



bedtools intersect -a $TIME_12 -b AB_geneIntersect.bed > AB_T12_variants.bam

bedtools intersect -a $TIME_12 -b TU_geneIntersect.bed > TU_T12_variants.bam

bedtools intersect -a $TIME_12 -b TL_geneIntersect.bed > TL_T12_variants.bam

bedtools intersect -a $TIME_12 -b WIK_geneIntersect.bed > WIK_T12_variants.bam


bedtools intersect -a $TIME_24 -b AB_geneIntersect.bed > AB_T24_variants.bam

bedtools intersect -a  $TIME_24 -b TU_geneIntersect.bed > TU_T24_variants.bam

bedtools intersect -a $TIME_24 -b TL_geneIntersect.bed > TL_T24_variants.bam

bedtools intersect -a $TIME_24 -b WIK_geneIntersect.bed > WIK_T24_variants.bam
 

#run featurecounts
featureCounts -a ${GTF} -o variantTranscriptsCounts.txt -T 24  --primary  -C  -p --countReadPairs AB_T0_variants.bam TU_T0_variants.bam TL_T0_variants.bam  WIK_T0_variants.bam AB_T6_variants.bam TU_T6_variants.bam TL_T6_variants.bam  WIK_T6_variants.bam AB_T12_variants.bam TU_T12_variants.bam TL_T12_variants.bam WIK_T12_variants.bam AB_T24_variants.bam  TU_T24_variants.bam TL_T24_variants.bam WIK_T24_variants.bam

