#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=512g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=abLine_ONT_PacBio


#Load the modules
module load verkko/2.1


#Run the reads extraction first

cd /data/okendojo/zebrafish/data/g3/assembly

#Run verkko
verkko -d  abLine_ONT_PacBio  --hifi /data/okendojo/zebrafish/data/g3/ab_inbreadPacBio/*.fastq.gz --nano /data/okendojo/zebrafish/data/g3/ontData/individual/*.fastq.gz --slurm  --lay-run 1 90 24
