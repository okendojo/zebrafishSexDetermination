#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=srf


#Load the modules
module load kmc
module load minimap2
module load k8

#Run the reads extraction first

cd  /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/pur

mkdir srf
cd srf

fish_6="../pur.fasta"

mkdir tmp_dir

kmc -fq -k151 -t16 -ci100 -cs100000 @input.txt count.kmc tmp_dir/

kmc_dump count.kmc count.txt

#assemble sat DNA
srf -p prefix count.txt > srf.fa


srfutils.js enlong srf.fa > srf.enlong.fa  # enlong short contigs for mapping

minimap2 -c -k 20 -N1000000 -f1000 -r100,100 -t 16 srf.enlong.fa ${fish_6} > fish6_srf.paf

srfutils.js paf2bed fish6_srf.paf > fish6_srf.bed  # generate non-redundant mapping regions

