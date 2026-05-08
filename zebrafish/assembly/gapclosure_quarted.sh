#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=GAP_closer


#Load modules
module load python/3.9 minimap2 trf cd-hit mummer/4.0.0beta2  R gnuplot/5.4.3 blast


cd  /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/alignments

quartet_gapfiller.py -d NHGRI_Fish6.fasta -t 24 -p gapclosing -g fish6.fasta  
