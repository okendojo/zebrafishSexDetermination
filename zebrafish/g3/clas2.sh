#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=90:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=tl_wik_classifyKmers

#activate conda
conda activate trio_binning

#load modules
module add gcc kmc
module load openmpi/4.1.3/gcc-11.3.0

cd /data/okendojo/zebrafish/data/g3/binning/tl_wik

classify-by-kmers \
    /data/okendojo/zebrafish/data/g3/pacBio/pacbio.fastq.gz \
    hapA_only_kmers.txt \
    hapB_only_kmers.txt \
    --haplotype-a-out-prefix pacbio_classified/maternal \
    --haplotype-b-out-prefix pacbio_classified/paternal \
    --unclassified-out-prefix pacbio_classified/unclassified


