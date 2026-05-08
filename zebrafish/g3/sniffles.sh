#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=210g
#SBATCH --ntasks-per-core=1
#SBATCH --time=10:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=sniffles

#load modules
module load sniffles

cd /data/okendojo/zebrafish/data/segmental_duplication/synteny/SV

#sniffles --input F6_Fish11.bam -t 24 --allow-overwrite --vcf F6F11.vcf

echo "processing AB...."

sniffles --input /data/okendojo/zebrafish/data/AB/polishing/syri/GRCz12tu_AB2.bam  -t 24 --allow-overwrite --vcf F6AB.vcf
