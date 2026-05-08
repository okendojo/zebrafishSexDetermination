#!/bin/sh

###############
##HEADER FOR BIOWULF SBATCH SUBMIT
##############
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --partition=norm
#SBATCH --mail-type=FAIL,END
#SBATCH --mem=100G
#SBATCH --gres=lscratch:400
#SBATCH --time=240:00:00
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=egapx
#######
## LOAD THE MODULES NEEDED TO RUN THE WORKFLOW
#######
module load python
module load nextflow
source /data/$USER/conda/etc/profile.d/conda.sh && source /data/$USER/conda/etc/profile.d/mamba.sh #source the conda env
conda activate egapx # activate conda env


cd /data/okendojo/zebrafish/data/annotation/fish6

python3 /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/eval/egapx/ui/egapx.py /data/okendojo/zebrafish/data/annotation/fish6/fish6.yaml -o NHGRI_ptx_Fish6 -e biowulf_cluster
