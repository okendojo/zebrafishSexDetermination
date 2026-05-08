#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=230g
#SBATCH --ntasks-per-core=1
#SBATCH --time=170:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=contigAlign


#Load the required modules
module load mafft
module load muscle

cd /data/okendojo/zebrafish/data/fish6/centromere/tmp

./FAMSA/famsa -t 32  -medoidtree  -v /data/okendojo/zebrafish/data/g3/rna_sequences/haplotype_gene_quantification/g3_variants/rna_variants/variant_calling/fourContigs.fasta  /data/okendojo/zebrafish/data/g3/rna_sequences/haplotype_gene_quantification/g3_variants/rna_variants/variant_calling/fourContigs_aln.fasta
