#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=90:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=strainspecifick-merinthecontigs

#load modules
module add meryl

cd /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls

#meryl-lookup -existence -sequence assembly.homopolymer-compressed.fasta -mers AB.only.meryl -output AB.txt

#meryl-lookup -existence -sequence assembly.homopolymer-compressed.fasta -mers TL.only.meryl -output TL.txt

#meryl-lookup -existence -sequence assembly.homopolymer-compressed.fasta -mers TU.only.meryl -output TU.txt

#meryl-lookup -existence -sequence assembly.homopolymer-compressed.fasta -mers WIK.only.meryl -output WIK.txt

meryl-lookup -existence -sequence /data/okendojo/zebrafish/data/g3/assembly/verkko_asm_run2/assembly.homopolymer-compressed.fasta -mers AB.only.meryl TU.only.meryl TL.only.meryl  WIK.only.meryl -output verkko_run1_contigKemrs.txt
