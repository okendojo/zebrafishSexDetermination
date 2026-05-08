#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=trf


#Load the modules
module load trf
module load kmc
module load minimap2
module load k8

#Run the reads extraction first

cd  /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/pur/srf/TRF-mod

#trf /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/pur/NHGRI_Fish6_cons.fasta  2 7 7 80 10 50 500 -l 10 -d -f -h  
 
./trf-mod  /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/pur/pur.fasta > PUR_trf.bed
