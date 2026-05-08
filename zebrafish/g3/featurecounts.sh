#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=200g
#SBATCH --ntasks-per-core=1
#SBATCH --time=10:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=featureCounts

#load modules
module load subread

cd /data/okendojo/zebrafish/data/g3/eQTL/fullTimePoints/preprocessing

GTF="/data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.111.gtf"

featureCounts -a ${GTF} -o grandparentsTranscriptCounts.txt -T 24  --primary  -C  -p --countReadPairs time_0/time_0.recal.bam time_6/time_6.recal.bam time_12/time_12.recal.bam time_24/time_24.recal.bam

#featureCounts -a ${GTF} -o variantTranscriptsCounts.txt -T 24  --primary  -C  -p --countReadPairs AB_T0_variants.bam TU_T0_variants.bam TL_T0_variants.bam  WIK_T0_variants.bam AB_T6_variants.bam TU_T6_variants.bam TL_T6_variants.bam  WIK_T6_variants.bam AB_T12_variants.bam TU_T12_variants.bam TL_T12_variants.bam WIK_T12_variants.bam AB_T24_variants.bam  TU_T24_variants.bam TL_T24_variants.bam WIK_T24_variants.bam 

