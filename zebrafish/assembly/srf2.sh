#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=srf2


#Load the modules
module load kmc
module load minimap2
module load k8

#Run the reads extraction first

cd /data/okendojo/zebrafish/data/segmental_duplication/srf


#fish_11="/data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/asm_polished/NHGRI_Fish6_cons.fasta"

fish_11="/data/okendojo/zebrafish/data/fish6/hifi/hifi.fasta"

kmc -fm -k151  -t16 -ci100  -cs100000 ${fish_11} fish11_count.kmc tmp_dir

kmc_dump fish11_count.kmc fish11_count.txt

srf -p prefix fish11_count.txt > fish11_srf.fa


srfutils.js enlong fish11_srf.fa > fish11_srf.enlong.fa  # enlong short contigs for mapping

minimap2 -c -k 20 -N1000000 -f1000 -r100,100 -t 16 fish11_srf.enlong.fa ${fish_11} > fish11_srf.paf

srfutils.js paf2bed fish11_srf.paf > fish11_srf.bed  # generate non-redundant mapping regions

