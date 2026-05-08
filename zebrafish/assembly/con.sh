#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=cons_fish6


#Load the modules
module load verkko/2.1
module load seqtk

#Run the reads extraction first

cd "/data/okendojo/zebrafish/data/fish6/asm"

#Run verkko

#seqtk subseq /data/okendojo/zebrafish/data/fish6/ontData/concatenated/f6_ont.fastq.gz "/data/okendojo/zebrafish/data/fish6/asm/f6_verkko_asm_v2.1/ontseqs.id" > "/data/okendojo/zebrafish/data/fish6/asm/f6_verkko_asm_v2.1/chr25filtered.fastq.gz"

verkko --assembly verkko_asm_f6 -d chrsTest25  --hifi /data/okendojo/zebrafish/data/fish6/hifi/m54313U_220817_024630.hifi_reads.fastq.gz --nano "/data/okendojo/zebrafish/data/fish6/ontData/concatenated/f6_ont.fastq.gz" --paths /data/okendojo/zebrafish/data/fish6/asm/consensusfile/consensus.gaf  --slurm

