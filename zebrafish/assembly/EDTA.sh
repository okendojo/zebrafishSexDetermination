#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=70:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=EDTA

#load modules
module load EDTA/2.2.2

cd /data/okendojo/zebrafish/data/AB/polishing/centrominer/edta

EDTA.pl --genome /data/okendojo/zebrafish/data/AB/polishing/syri/AB.fa  --cds /data/okendojo/zebrafish/data/AB/polishing/centrominer/cds.faa --sensitive 1 --overwrite 1 --anno 1 --threads 24 --force 1
