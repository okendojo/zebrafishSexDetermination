#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=teloExplorer


cd /data/okendojo/zebrafish/data/fish6/centromere

quartet_teloexplorer.py -i  fish6.fasta -p fish_6Telo 

