#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=stringtie

#load modules
module load stringtie
module load bedtools

cd /data/okendojo/zebrafish/data/annotation/fish6/others

stringtie Fish6.pri.bam -o GRCz112tu_transcripts.gtf -v -p 24

bedtools genomecov -ibam  Fish6.pri.bam -bg > GRCz12tu_coverage.bedgraph


