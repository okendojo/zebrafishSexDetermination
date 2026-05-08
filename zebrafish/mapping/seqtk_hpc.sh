#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=seqtk_hpc

#load the modules
module load seqtk

cd /data/okendojo/zebrafish/data/AB/batches_ont

seqtk seq -A ont.fastq.gz > ont.fasta

seqtk hpc ont.fasta > ont_hpc.fasta

rm ont.fasta 

