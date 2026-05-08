#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=F2_RNA_meryl_lookup

# Load meryl module
module add meryl

cd /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls

meryl-lookup -existence -sequence  /data/okendojo/zebrafish/data/g3/assembly/asm_g3fullont/assembly.homopolymer-compressed.fasta -mers AB.only.meryl TU.only.meryl  TL.only.meryl  WIK.only.meryl -output f2_contigKemrs.txt

#meryl-lookup -bed -sequence /data/okendojo/zebrafish/data/g3/ontData/combined/ont.fasta  -mers AB.only.meryl/ TU.only.meryl/ TL.only.meryl/ WIK.only.meryl/ -output ont_readsMatch.bed  -labels AB TU TL WIK 


