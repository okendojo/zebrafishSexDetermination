#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=64g
#SBATCH --ntasks-per-core=1
#SBATCH --time=96:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov

#load the module
module load bwa
module load samtools
module load dragmap/1.2.1
module load GATK
module load picard

cd /data/okendojo/paradisfishProject/snv/dragmapping


#Add read groups else haplotype caller will not run
#java -jar $PICARDJARPATH/picard.jar AddOrReplaceReadGroups -I macOpe2.bam  -O macOpe2_RG.bam -RGID 4 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM 20 -SO coordinate --CREATE_INDEX true


#mark duplicates
#java -jar $PICARDJARPATH/picard.jar  MarkDuplicates -I macOpe2_RG.bam -O macOpe2_RG.dedup.bam -M marked_dup_metrics.txt


#call variants
gatk  HaplotypeCaller  \
   -R ../assembly/macOpe2Assembly.fasta \
   -I macOpe2.sorted.bam \
   -O mc_output.g.vcf.gz \
   -ERC GVCF 






