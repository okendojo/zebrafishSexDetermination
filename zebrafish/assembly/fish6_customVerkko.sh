#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=fish6_verkko


#Load the modules
module load verkko/2.1

#Run the reads extraction first

cd /data/okendojo/zebrafish/data/fish6/asm


#Run verkko
verkko -d  asm_f6_Herro_newTest  --hifi /data/okendojo/zebrafish/data/fish6/hifi/m54313U_220817_024630.hifi_reads.fastq.gz  --nano /data/okendojo/zebrafish/data/fish6/ontData/concatenated/*.fa --slurm

