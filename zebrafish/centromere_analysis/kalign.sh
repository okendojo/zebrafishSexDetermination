#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=170:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=kalign


#Load the required modules
module load kalign


cd /data/okendojo/zebrafish/data/g3/rna_sequences/haplotype_gene_quantification

fasta="/data/okendojo/zebrafish/data/g3/rna_sequences/haplotype_gene_quantification/g3_variants/rna_variants/variant_calling/hapTagsrenamed.fasta"

kalign  -i ${fasta} -f fasta --type dna --nthreads 24 > hapTagsrenamed_aln.fasta
