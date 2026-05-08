#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=230g
#SBATCH --ntasks-per-core=1
#SBATCH --time=140:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=G3_verkkoABTU


#Load the modules

module load verkko

cd /data/okendojo/zebrafish/data/g3/assembly

verkko -d asm_abtu_2 --grid --hifi /data/okendojo/zebrafish/data/g3/output/meryl/bins/abtu/hifi/ABTU_PacBio.fasta.gz  --nano /data/okendojo/zebrafish/data/g3/output/meryl/bins/abtu/ont/ABTU_ONT.fasta.gz --hap-kmers /data/okendojo/zebrafish/data/g3/output/blobplots/verkko_gpcompressedmeryls/abtu_child/TU.k30.hapmer.meryl /data/okendojo/zebrafish/data/g3/output/blobplots/verkko_gpcompressedmeryls/abtu_child/AB.k30.hapmer.meryl trio  
