#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=64g
#SBATCH --ntasks-per-core=1
#SBATCH --time=96:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=eQTL

#load the module
module load bwa
module load samtools
module load GATK
module load picard
module add bcftools

cd /data/okendojo/zebrafish/data/g3/rna_sequences/gene_quantification/eQTL 

#gatk SelectVariants  -R /data/okendojo/zebrafish/data/g3/rna_sequences/gene_quantification/rnaseq_result/star_salmon/Danio_rerio.GRCz11.dna.primary_assembly.fa  --variant /data/okendojo/zebrafish/data/g3/rna_sequences/gene_quantification/rnaVariants.vcf --restrict-alleles-to BIALLELIC -select 'vc.getHetCount()==1' --select-type-to-include SNP -O dna_variants.selected.vcf

#bcftools norm --rm-dup all dna_variants.selected.vcf | bgzip > out.vcf.gz

gatk ASEReadCounter --input /data/okendojo/zebrafish/data/g3/rna_sequences/gene_quantification/rnaseq_result/star_salmon/time_0_1.markdup.sorted.bam -I	/data/okendojo/zebrafish/data/g3/rna_sequences/gene_quantification/rnaseq_result/star_salmon/time_6_1.markdup.sorted.bam  -I /data/okendojo/zebrafish/data/g3/rna_sequences/gene_quantification/rnaseq_result/star_salmon/time_12_1.markdup.sorted.bam -I /data/okendojo/zebrafish/data/g3/rna_sequences/gene_quantification/rnaseq_result/star_salmon/time_24_1.markdup.sorted.bam --variant out.vcf.gz --output eqtlResults.table --min-base-quality 10 --min-mapping-quality 10 --reference /data/okendojo/zebrafish/data/g3/rna_sequences/gene_quantification/rnaseq_result/star_salmon/Danio_rerio.GRCz11.dna.primary_assembly.fa --sequence-dictionary /data/okendojo/zebrafish/data/g3/rna_sequences/gene_quantification/rnaseq_result/star_salmon/Danio_rerio.GRCz11.dna.primary_assembly.dict 
