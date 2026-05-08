#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=36:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=quast

#Activate mamba environment
module load augustus
module add quast

cd /data/okendojo/zebrafish/data/g3/assembly

WIK="/data/okendojo/zebrafish/data/g3/assembly/wik_asm/assembly.fasta"
TL="/data/okendojo/zebrafish/data/g3/assembly/tl_asm/assembly.fasta"
AB="/data/okendojo/zebrafish/data/g3/assembly/ab_asm/assembly.fasta"
TU="/data/okendojo/zebrafish/data/g3/assembly/tu_asm/assembly.fasta"


quast.py -o quast_gene -r /data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.dna.primary_assembly.fa --est-ref-size 1373471384 --debug --fragmented  -l "WIK, TL, AB, TU" --space-efficient --eukaryote -t 30 --min-identity 95.0 --plots-format pdf --glimmer --circos  $WIK $TL $AB $TU   

