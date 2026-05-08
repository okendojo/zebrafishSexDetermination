#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=360g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=cellrangercount

#load modules
#ml cellranger
ml cellranger-arc

cd /data/okendojo/zebrafish/data/g3/scanalysis

echo "Now working on GRCz12tu2......"

cellranger-arc count --id=GRCz12tu_arc_noformatasm --reference=/data/okendojo/zebrafish/data/g3/scanalysis/genomeasm/mkref/GRCz12tu --libraries=/data/okendojo/zebrafish/data/g3/scanalysis/library/sample_multiome_csv.csv --jobmode slurm --description GRCz12tu

#echo "Now running GRCz11 assembly..."

#cellranger-arc count --id=GRCz11_arc_force --reference=/data/okendojo/zebrafish/data/g3/scanalysis/ref/grcz11/GRCz11 --libraries=/data/okendojo/zebrafish/data/g3/scanalysis/library/sample_multiome_csv.csv --jobmode slurm --description GRCz11  --force-cells=10000
