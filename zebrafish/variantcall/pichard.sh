#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=64g
#SBATCH --ntasks-per-core=1
#SBATCH --time=96:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov

module load picard

cd /data/okendojo/paradisfishProject/snv

java -jar $PICARDJARPATH/picard.jar AddOrReplaceReadGroups -I macOpe2.bam  -O macOpe2_RG.bam -RGID 4 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM 20 -SO coordinate --CREATE_INDEX true

