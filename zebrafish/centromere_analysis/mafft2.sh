#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=200g
#SBATCH --ntasks-per-core=1
#SBATCH --time=170:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=mafft_alignment


#Load the required modules
module load mafft
module load muscle

cd /data/okendojo/zebrafish/data/g3/rna_sequences/haplotype_gene_quantification

fasta="/data/okendojo/zebrafish/data/g3/rna_sequences/haplotype_gene_quantification/g3_variants/rna_variants/variant_calling/fourContigs.fasta"


mafft --thread 24 --maxiterate 100 ${fasta} > fourcontigs_mafft_aln.fasta  
