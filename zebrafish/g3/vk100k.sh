#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=100kG3_verkko


#Load the modules
module load seqtk
module load verkko/2.2


#Run the reads extraction first

cd  /data/okendojo/zebrafish/data/g3/assembly


#Run verkko
verkko -d asm_g3_100k  --hifi /data/okendojo/zebrafish/data/g3/pacBio/*.fastq.gz --nano  /data/okendojo/zebrafish/data/g3/ontData/combined/filt100k.fasta --slurm
