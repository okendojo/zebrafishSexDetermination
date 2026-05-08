#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=80g
#SBATCH --ntasks-per-core=1
#SBATCH --time=96:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#move to the working directory

#Load the modules
module load merqury/1.3

cd /data/okendojo/zebrafish/data/g3/hsp 

meryl count k=21 output hifi.meryl /data/okendojo/zebrafish/data/g3/pacBio/*.fastq.gz 

