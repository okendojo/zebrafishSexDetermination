#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=64g
#SBATCH --ntasks-per-core=1
#SBATCH --time=96:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov

# Workflow is adopted from: https://gatk.broadinstitute.org/hc/en-us/articles/4407897446939--How-to-Run-germline-single-sample-short-variant-discovery-in-DRAGEN-mode#article-comments
#Load the required module
module load GATK
module load dragmap
#module load singularity
#move to the right dir
cd  /data/okendojo/paradisfishProject/snv


#Build the DRAGEN STR Model from the Aligned Reads

gatk CalibrateDragstrModel -R assembly/macOpe2Assembly.fasta -I macOpe2_RG.bam -str str_table.tsv -O dragstr_model.txt

#Variant Calling and Filtering
#gatk HaplotypeCaller -R ${REF} -I macOpe2.sam -O macOpe2_snvs.vcf --dragen-mode true --dragstr-params-path dragstr_model.txt

#Hard Filter Variants
#gatk VariantFiltration -V macOpe2_snvs.vcf --filter-expression "QUAL < 10.4139" --filter-name "DRAGENHardQUAL" -O macOpe2_snv_filtered.vcf
















