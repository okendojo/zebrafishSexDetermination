#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=ab_binning

#load modules
module add gcc kmc python
module load openmpi/4.1.3/gcc-11.3.0

#source conda

#cd /data/okendojo/zebrafish/data/g3/pacBio

#cat *.fastq.gz > hifi.fastq.gz

source /data/$USER/conda/etc/profile.d/conda.sh && source /data/$USER/conda/etc/profile.d/mamba.sh

#activate conda
conda activate trio_binning

cd /data/okendojo/zebrafish/data/ab_asm/binning/dfraw

classify-by-kmers \
    /data/okendojo/zebrafish/data/g3/pacBio/hifi.fastq.gz \
    /data/okendojo/zebrafish/data/ab_asm/binning/finduniqkmers/hapA_only_kmers.txt \
    /data/okendojo/zebrafish/data/ab_asm/binning/finduniqkmers/hapB_only_kmers.txt \
    --haplotype-a-out-prefix classified/AB_mat \
    --haplotype-b-out-prefix classified/TU_pat \
    --unclassified-out-prefix classified/unclassified


