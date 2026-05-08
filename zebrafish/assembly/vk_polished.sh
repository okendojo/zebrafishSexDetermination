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

#Run verkko

cd  /data/okendojo/zebrafish/data/AB/asm

#Run verkko
verkko -d 2_curated_ont_ul_asm --paths /data/okendojo/zebrafish/data/AB/asm/ont_ul_asm/asm_graph_map/2_remaining.gaf  --assembly ont_ul_asm --hifi /data/okendojo/zebrafish/data/AB/hifi/*.gz --nano /data/okendojo/zebrafish/data/AB/batches_ont/ont_ul.fastq.gz --grid --lay-run 16 100 96 --cns-run 16 200 36
