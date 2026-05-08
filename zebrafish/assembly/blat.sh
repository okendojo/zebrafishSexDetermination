#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=blat

module load blat

cd /data/okendojo/zebrafish/data/fish6/centromere/teltools 

blat ../fish6.fasta  telomeric_specific_reads.fasta -t=DNA -q=DNA -minScore=50 -minIdentity=95  f6blat_output.psl 
