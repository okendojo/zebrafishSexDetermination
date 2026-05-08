#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=Verkko


#Load the modules
module load verkko/2.1
module load seqkit

cd /data/okendojo/zebrafish/data/g3/assembly

#Run verkko
verkko -d vk_ONT_HiFideduped  --hifi /data/okendojo/zebrafish/data/g3/pacBio/*.fastq.gz --nano /data/okendojo/zebrafish/data/g3/ontData/combined/ont.fastq.gz --slurm
