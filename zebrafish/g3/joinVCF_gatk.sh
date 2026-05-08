#!/bin/bash
#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=200g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=transcriptvar

#load modules
module add GATK
module add samtools

# Define paths and filenames
REFERENCE_GENOME="/data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.dna.primary_assembly.fa"
OUTPUT_DIR="/data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/bamFiles/GATKCall"
BAM_DIR="/data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/bamFiles"

cd $BAM_DIR

for bam in *_sorted.bam ; do samtools index $bam ; done 

# Run individual bam files
gatk HaplotypeCaller -R $REFERENCE_GENOME -I $BAM_DIR/AB_T0_variants.bam_sorted.bam -O $OUTPUT_DIR/AB_T0_variants.vcf

gatk HaplotypeCaller -R $REFERENCE_GENOME -I $BAM_DIR/AB_T12_variants.bam_sorted.bam -O $OUTPUT_DIR/AB_T12_variants.vcf

gatk HaplotypeCaller -R $REFERENCE_GENOME -I $BAM_DIR/AB_T24_variants.bam_sorted.bam -O $OUTPUT_DIR/AB_T24_variants.vcf

gatk HaplotypeCaller -R $REFERENCE_GENOME -I $BAM_DIR/AB_T6_variants.bam_sorted.bam -O $OUTPUT_DIR/AB_T6_variants.vcf

gatk HaplotypeCaller -R $REFERENCE_GENOME -I $BAM_DIR/TL_T0_variants.bam_sorted.bam -O $OUTPUT_DIR/TL_T0_variants.vcf

gatk HaplotypeCaller -R $REFERENCE_GENOME -I $BAM_DIR/TL_T12_variants.bam_sorted.bam -O $OUTPUT_DIR/TL_T12_variants.vcf

gatk HaplotypeCaller -R $REFERENCE_GENOME -I $BAM_DIR/TL_T24_variants.bam_sorted.bam -O $OUTPUT_DIR/TL_T24_variants.vcf

gatk HaplotypeCaller -R $REFERENCE_GENOME -I $BAM_DIR/TL_T6_variants.bam_sorted.bam -O $OUTPUT_DIR/TL_T6_variants.vcf

gatk HaplotypeCaller -R $REFERENCE_GENOME -I $BAM_DIR/TU_T0_variants.bam_sorted.bam -O $OUTPUT_DIR/TU_T0_variants.vcf

gatk HaplotypeCaller -R $REFERENCE_GENOME -I $BAM_DIR/TU_T12_variants.bam_sorted.bam -O $OUTPUT_DIR/TU_T12_variants.vcf

gatk HaplotypeCaller -R $REFERENCE_GENOME -I $BAM_DIR/TU_T24variants.bam_sorted.bam -O $OUTPUT_DIR/TU_T24_variants.vcf

gatk HaplotypeCaller -R $REFERENCE_GENOME -I $BAM_DIR/TU_T6_variants.bam_sorted.bam -O $OUTPUT_DIR/TU_T6_variants.vcf

gatk HaplotypeCaller -R $REFERENCE_GENOME -I $BAM_DIR/WIK_T0_variants.bam_sorted.bam -O $OUTPUT_DIR/WIK_T0_variants.vcf

gatk HaplotypeCaller -R $REFERENCE_GENOME -I $BAM_DIR/WIK_T12_variants.bam_sorted.bam -O $OUTPUT_DIR/WIK_T12_variants.vcf

gatk HaplotypeCaller -R $REFERENCE_GENOME -I $BAM_DIR/WIK_T24_variants.bam_sorted.bam -O $OUTPUT_DIR/WIK_T24_variants.vcf

gatk HaplotypeCaller -R $REFERENCE_GENOME -I $BAM_DIR/WIK_T6_variants.bam_sorted.bam -O $OUTPUT_DIR/WIK_T6_variants.vcf



