#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov

#load the module
module load VEP


cd /data/okendojo/zebrafish/data/g3/rna_sequences/haplotype_gene_quantification/rna_variants/variant_calling/variants 


for vcf in *.vcf.gz ; do vep -i $vcf -o ${vcf}_annot.out --offline --cache --force_overwrite --dir ${VEP_CACHEDIR} --species zebrafish --fasta $VEP_CACHEDIR/zebrafish.fa ; echo "processing file $vcf:" ; done 

