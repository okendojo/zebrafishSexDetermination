#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=230g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=G3asm_meryl_lookup

# Load meryl module
module add meryl

cd /data/okendojo/zebrafish/data/g3/illumina/meryldbs

meryl-lookup -existence -sequence /data/okendojo/zebrafish/data/g3/ontData/combined/ont.fasta  -mers AB.only.meryl/ TU.only.meryl/ TL.only.meryl/ WIK.only.meryl/ -output ont_fasta.bed  #-labels AB TU TL WIK

