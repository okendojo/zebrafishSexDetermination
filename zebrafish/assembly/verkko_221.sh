#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=230g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=AB_verkko


#Load the modules
module load verkko
module load python/3.10

#Run the reads extraction first

cd  /data/okendojo/zebrafish/data/AB/asm

ref="/data/ShawnBurgessLab/javan_duncan_multiome/grcz12tuannotations/GRCz12tu.fasta"

#Run verkko
verkko -d ont_ul_asm --hifi /data/okendojo/zebrafish/data/AB/hifi/*.gz --nano /data/okendojo/zebrafish/data/AB/batches_ont/ont_ul.fastq.gz --ref ${ref} --grid --lay-run 16 100 96 --cns-run 16 200 36
