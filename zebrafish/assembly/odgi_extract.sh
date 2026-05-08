#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=750g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH	--gres=lscratch:1000
#SBATCH --job-name=extract

cd /data/okendojo/zebrafish/data/g3/sex_project/danio_pg

#odgi extract -i  danio_pg.full.og -o tmp_mhc.og -b gene_mhc.bed -c 0 -E --threads 24 -P

#odgi paths -i tmp_mhc.og -L -t 20 -P

#odgi stats -i tmp_mhc.og -S -t 20 -P

odgi view -i tmp_mhc.og -g -t 20 -P > chr25.gfa
