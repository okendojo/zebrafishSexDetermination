#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=750g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --gres=lscratch:1000
#SBATCH --job-name=odgi_viz

cd /data/okendojo/zebrafish/data/g3/sex_project/assemblies/chr19/Danio_pg1

odgi procbed -i tmp.gfa -b mhc.bed > adj.bed 
