#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=170:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=SRF_F11


#Load the required modules
module add k8
module add kmc
module add minimap2

DIR="/data/okendojo/zebrafish/data/centromere_analysis/satelliteRepFinder/fish11"

#change to the working directory
cd ${DIR}

mkdir -p tmp_dir 

#count the kmers using KMC
kmc -fq -k171 -t16 -ci100 -m200 -cs100000  -sf24 -sp24 -sr24 @hifi.txt count.kmc tmp_dir

kmc_dump count.kmc count.txt

#run SRF to get contigs
srf -p Fish11 count.txt > srf.fa 

srfutils.js enlong srf.fa > srf.enlong.fa  # enlong short contigs for mapping

minimap2 -c -N1000000 -f1000 -r100,100 -t16 srf.enlong.fa /data/Zebrafish_T2T/fish11/resolved_assembly/fish11_merfin_polished.fasta > srf.paf

srfutils.js paf2bed srf.paf > srf.bed  # generate non-redundant mapping regions

srfutils.js bed2abun srf.bed > srf_abundance.len # calculate abundance of each srf contig
