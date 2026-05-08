#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=f6_telogator


#Load the modules
module load minimap2
module load winnowmap
module load python/3.9


cd /data/okendojo/zebrafish/data/fish6/centromere

python telogator2/telogator2.py -i /data/okendojo/zebrafish/data/fish6/centromere/fish6.fasta -o f6_telogator -p 24 --muscle bin/muscle -r hifi -tt 0.400 -ts 0.250 -n 4 
