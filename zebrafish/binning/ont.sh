#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=ont_binning

#load modules
module add python gcc kmc python
module load openmpi/4.1.3/gcc-11.3.0

#source conda
source /data/$USER/conda/etc/profile.d/conda.sh && source /data/$USER/conda/etc/profile.d/mamba.sh

#activate conda
conda activate trio_binning

cd /data/okendojo/zebrafish/data/ab_asm/binning/dfraw

classify-by-kmers \
    /data/okendojo/zebrafish/data/g3/ontData/combined/comb_ont.fasta \
    /data/okendojo/zebrafish/data/ab_asm/binning/finduniqkmers/hapA_only_kmers.txt \
    /data/okendojo/zebrafish/data/ab_asm/binning/finduniqkmers/hapB_only_kmers.txt \
    --haplotype-a-out-prefix  ont/AB_mat \
    --haplotype-b-out-prefix  ont/TU_pat \
    --unclassified-out-prefix  ont/unclassified


