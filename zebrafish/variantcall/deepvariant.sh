#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=64g
#SBATCH --ntasks-per-core=1
#SBATCH --time=96:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=deepvariant

module load samtools

cd /data/okendojo/zebrafish/data/g3/rna_sequences/gene_quantification/rnaseq_result/star_salmon

#samtools view -S -b  alignment.sam > illumeAlignment.bam
#samtools sort -o illumeAlignment_sorted.bam illumeAlignment.bam 
#samtools index illumeAlignment_sorted.bam 

samtools merge -b bamList.fofn -O BAM mergedBamFile.bam -f
samtools sort -o mergedBamFile_sorted.bam mergedBamFile.bam

samtools index mergedBamFile_sorted.bam 

rm mergedBamFile.bam 

./_submit_deepvariant.sh Danio_rerio.GRCz11.dna.primary_assembly.fa  mergedBamFile_sorted.bam  WGS RNA_variants
