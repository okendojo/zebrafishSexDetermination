#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=stampyAlignment


#move to the directory containing the data
cd /data/okendojo/zebrafish/data/G3/TU_AB_alignment

#Set the path variables for the genome index and the harsh table
index=/data/okendojo/zebrafish/refGenome/index/genome
hashtable=/data/okendojo/zebrafish/refGenome/index/genome

TU_BAM=/data/okendojo/zebrafish/data/G3/TU_AB_alignment/TU.bam	
AB_BAM=/data/okendojo/zebrafish/data/G3/TU_AB_alignment/AB.bam


#Run sampy to select good reads only
stampy.py -g ${index} -h ${hashtable} -t 32 --bamkeepgoodreads -f sam -o TU_stampyCleaned.sam -M ${TU_BAM} 

stampy.py -g ${index} -h ${hashtable} -t 32 --bamkeepgoodreads -f sam -o AB_stampCleaned.sam -M ${AB_BAM}  
