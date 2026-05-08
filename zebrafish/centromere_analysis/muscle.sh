#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=170:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=centromere_analysis


#Load the required modules
module load mafft
module load muscle

cd /data/okendojo/zebrafish/data/fish6/centromere/tmp

fasta="/data/okendojo/zebrafish/data/g3/rna_sequences/haplotype_gene_quantification/g3_variants/rna_variants/variant_calling/hapTagsrenamed.fasta"

muscle -align ${fasta} -output topContigs_aln.afa
