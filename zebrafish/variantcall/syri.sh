#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=64g
#SBATCH --ntasks-per-core=1
#SBATCH --time=96:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov

#Source conda & mamba
# Workflow adopted from: https://github.com/schneebergerlab/plotsr#visualising-tracks

source /data/$USER/conda/etc/profile.d/conda.sh && source /data/$USER/conda/etc/profile.d/mamba.sh

# activate conda env
mamba activate plotsr

cd /data/okendojo/paradisfishProject/genoComparison/

syri -c macOpe2.sorted.bam -r betta2.fasta -q macope.fasta -F B --prefix Betta_macOpe2 --dir syri  
