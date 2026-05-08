#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=360g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=WIK


#Load the modules
module load verkko


#Run the reads extraction first

cd  /data/okendojo/zebrafish/data/g3/assembly

#Run verkko
verkko -d wik_asm100k  --hifi /data/okendojo/zebrafish/data/ab_asm/uncompressedGPmeryls/uncompressed/fasta_ont_hifi/wik_hifi.fasta --nano /data/okendojo/zebrafish/data/ab_asm/uncompressedGPmeryls/uncompressed/fasta_ont_hifi/wik_ont100k.fasta  --slurm
