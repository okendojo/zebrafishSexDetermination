#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=hicat


#Load the modules
module load python

cd /data/okendojo/zebrafish/data/fish6/centromere/HiCAT

python HiCAT.py -i ../fish11_centromera.fasta -t testdata/AlphaSat.fa -th 20 -sn 20 -o fish11

python HiCAT.py -i ../fish6.fasta -t testdata/AlphaSat.fa -th 20 -sn 20 -o fish6
