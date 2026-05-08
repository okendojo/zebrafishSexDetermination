#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=edta


#Load the modules
module load genometools
module load repeatmasker repeatmodeler

cd /data/okendojo/zebrafish/data/fish6/asm/t2t/repeat_analysis

EDTA.pl --genome /data/okendojo/zebrafish/data/fish6/asm/t2t/NHGRI_Fish6.fasta --cds zfish.fasta --threads 24 --anno 1
