#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=170:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=famsa_newick


#Load the required modules
module load mafft
module load muscle

cd /data/okendojo/zebrafish/data/fish6/centromere/tmp

./FAMSA/famsa -t 24  -gt_export -v -medoidtree f11.fasta  fish11_fish6_cent_aln.newick 
