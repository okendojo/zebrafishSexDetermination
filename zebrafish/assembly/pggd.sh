#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=750g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH	--gres=lscratch:1000
#SBATCH --job-name=pggd_alignment


#Load the modules
ml pggb

cd /data/okendojo/zebrafish/data/g3/sex_project/pggd

pggb -i /vf/users/okendojo/zebrafish/data/g3/sex_project/assemblies/asm.fasta -o output -n 12 -Z -K 21 -m -t 16 -p 90 -s 5k
