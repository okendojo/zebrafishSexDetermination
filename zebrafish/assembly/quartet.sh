#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=60:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=centrominer

#Load modules
module load python/3.9 minimap2 trf cd-hit mummer/4.0.0beta2  R gnuplot/5.4.3 blast


cd /data/okendojo/zebrafish/data/AB/polishing/centrominer

quartet_centrominer.py -i /data/okendojo/zebrafish/data/AB/polishing/syri/AB.fa --gene /data/ShawnBurgessLab/javan_duncan_multiome/grcz12tuannotations/GRCz12tu_genomic.gtf --TE /data/okendojo/zebrafish/data/segmental_duplication/edta/NHGRI_Fish6_cons.fasta.mod.EDTA.TEanno.gff3 -t 24 -p GRCz12ab -s 0.95 -l 300000  
