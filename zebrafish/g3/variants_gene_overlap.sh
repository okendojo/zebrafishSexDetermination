#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=timePoints_var_overlap

cd  /data/okendojo/zebrafish/data/g3/rna_sequences/timepointsQuant/f0_genes 

./prcs.sh
