#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=750g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH	--gres=lscratch:1000
#SBATCH --job-name=mhc_pggd_alignment


#Load the modules
ml pggb

cd /data/okendojo/zebrafish/data/g3/sex_project/assemblies/chr19

pggb -i asm.fasta -o mhc1uka_pggd -n 12 -Z -K 21 -m -t 16 -p 90 -s 5k
