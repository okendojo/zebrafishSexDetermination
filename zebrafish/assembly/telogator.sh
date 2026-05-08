#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=telogatorr


#Load the modules
module load minimap2
module load winnowmap
module load python/3.9


cd /data/okendojo/zebrafish/data/fish6/centromere

python telogator2/telogator2.py -i /data/okendojo/zebrafish/data/fish6/hifi/m54313U_220817_024630.hifi_reads.fastq.gz -o tgator_telomeres -p 24 --muscle bin/muscle -r hifi -tt 0.400 -ts 0.250 -n 4 
