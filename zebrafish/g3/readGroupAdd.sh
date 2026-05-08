#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=64g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=addReadGroup

#load the module
module load bwa
module load samtools
module load GATK
module load picard
module load bamtools

#move to the working directory

cd /data/okendojo/zebrafish/data/g3/eQTL/BAMs

samtools addreplacerg -r '@RG\tID:samplename\tSM:samplename' 22CMMHLT3_19570453_S61_L003/Aligned.sortedByCoord.out.bam -o 22CMMHLT3_19570453_S61_L003.bam

samtools addreplacerg -r '@RG\tID:samplename\tSM:samplename' 22CMMHLT3_19570455_S62_L003/Aligned.sortedByCoord.out.bam -o 22CMMHLT3_19570455_S62_L003.bam 
	
samtools addreplacerg -r '@RG\tID:samplename\tSM:samplename' 22CMMHLT3_19570457_S63_L003/Aligned.sortedByCoord.out.bam -o 22CMMHLT3_19570457_S63_L003.bam 

samtools addreplacerg -r '@RG\tID:samplename\tSM:samplename' 22CMMHLT3_19570459_S64_L003/Aligned.sortedByCoord.out.bam -o 22CMMHLT3_19570459_S64_L003.bam

samtools addreplacerg -r '@RG\tID:samplename\tSM:samplename' 22CMMHLT3_19570461_S65_L003/Aligned.sortedByCoord.out.bam -o 22CMMHLT3_19570461_S65_L003.bam

samtools addreplacerg -r '@RG\tID:samplename\tSM:samplename' 22CMMHLT3_19570463_S66_L003/Aligned.sortedByCoord.out.bam -o 22CMMHLT3_19570463_S66_L003.bam

samtools addreplacerg -r '@RG\tID:samplename\tSM:samplename' 22CMMHLT3_19570465_S67_L003/Aligned.sortedByCoord.out.bam -o 22CMMHLT3_19570465_S67_L003.bam

samtools addreplacerg -r '@RG\tID:samplename\tSM:samplename' 22CMMHLT3_19570467_S68_L003/Aligned.sortedByCoord.out.bam -o 22CMMHLT3_19570467_S68_L003.bam

samtools addreplacerg -r '@RG\tID:samplename\tSM:samplename' 22CMMHLT3_19570469_S69_L003/Aligned.sortedByCoord.out.bam -o 22CMMHLT3_19570469_S69_L003.bam

samtools addreplacerg -r '@RG\tID:samplename\tSM:samplename' 22CMMHLT3_19570471_S70_L003/Aligned.sortedByCoord.out.bam -o 22CMMHLT3_19570471_S70_L003.bam

for bam in *.bam; do
	
	samtools index $bam

done
