#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=200g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=rmDUP

#load modules
module add nextflow
module add singularity
module load picard


cd /data/okendojo/zebrafish/data/g3/eQTL/fullTimePoints/preprocessing

java -Xmx24g -XX:ParallelGCThreads=16 -jar $PICARDJARPATH/picard.jar MarkDuplicates -I /data/okendojo/zebrafish/data/g3/eQTL/fullTimePoints/preprocessing/f2_rnalignments.bam -M dupMat.txt -O f2_rnalignments_dedup.bam --CREATE_INDEX true --REMOVE_DUPLICATES true 
