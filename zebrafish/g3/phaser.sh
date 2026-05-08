#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=200g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=ASEP

#load modules
module load phaser

cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants

BAM_T0="/data/okendojo/zebrafish/data/g3/eQTL/fullTimePoints/preprocessing/time_0/time_0.recal.bam"
BAM_T6="/data/okendojo/zebrafish/data/g3/eQTL/fullTimePoints/preprocessing/time_6/time_6.recal.bam"
BAM_T12="/data/okendojo/zebrafish/data/g3/eQTL/fullTimePoints/preprocessing/time_12/time_12.recal.bam"
BAM_T24="/data/okendojo/zebrafish/data/g3/eQTL/fullTimePoints/preprocessing/time_24/time_24.recal.bam"


phaser.py --vcf f0_annotatedVariants.vcf.gz --bam ${BAM_T0} --sample AB_AB_strain --mapq 255 --threads 24 --o AB_T0_hapCount.txt --baseq 10 --write_vcf 1 --paired_end 1

phaser.py --vcf f0_annotatedVariants.vcf.gz --bam ${BAM_T0} --sample TU_TU_strain --mapq 255 --threads 24 --o TU_T0_hapCount.txt --baseq 10 --write_vcf 1 --paired_end 1

phaser.py --vcf f0_annotatedVariants.vcf.gz --bam ${BAM_T0} --sample TL_TL_strain --mapq 255 --threads 24 --o TL_T0_hapCount.txt --baseq 10 --write_vcf 1 --paired_end 1

phaser.py --vcf f0_annotatedVariants.vcf.gz --bam ${BAM_T0} --sample WIK_WIK_strain --mapq 255 --threads 24 --o WIK_T0_hapCount.txt --baseq 10 --write_vcf 1 --paired_end 1



