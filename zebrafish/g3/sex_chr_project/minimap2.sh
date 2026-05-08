#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=minimap2

ml minimap2

cd /data/okendojo/zebrafish/data/g3/sex_project/scinkd_analysis/WW
minimap2 -x asm5 -t 24 -c --eqx --secondary=no WW.hap1.fasta.gz WW.hap2.fasta.gz -o WW.paf 

cd /data/okendojo/zebrafish/data/g3/sex_project/scinkd_analysis/ZW

minimap2 -x asm5 -t 24 -c --eqx --secondary=no ZW.hap1.fasta.gz ZW.hap2.fasta.gz -o ZW.paf 

cd /data/okendojo/zebrafish/data/g3/sex_project/scinkd_analysis/ZZ

bgzip ZZ.hap1.fasta  
bgzip ZZ.hap2.fasta

minimap2 -x asm5 -t 24 -c --eqx --secondary=no ZZ.hap1.fasta.gz  ZZ.hap2.fasta.gz -o ZZ.paf 
