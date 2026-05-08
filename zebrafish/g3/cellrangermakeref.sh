#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=360g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=makeref

#load modules
#ml cellranger
ml cellranger-arc

cd  /data/okendojo/zebrafish/data/g3/scanalysis/genomeasm/mkref

cellranger-arc mkref --config=/data/okendojo/zebrafish/data/g3/scanalysis/ref/GRCz12tu.config --nthreads=24 --memgb=100

#cellranger mkref --genome GRCz11 --fasta /data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.dna.primary_assembly.fa --jobmode slurm --genes /data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.111.gtf --nthreads 16 --output-dir GRCz11


#cd /data/okendojo/zebrafish/data/g3/scanalysis/ref/grcz11

#cellranger-arc mkref --config=/data/okendojo/zebrafish/data/g3/scanalysis/ref/GRCz11.config --nthreads=24 --memgb=100

#cellranger mkref --genome GRCz12 --fasta /data/okendojo/zebrafish/data/g3/sex_project/assemblies/NHGRI_Fish6_cons.fasta --jobmode slurm --genes /data/okendojo/zebrafish/data/annotation/grcz12tu_lessCurated.gtf  --nthreads 16 --output-dir GRCz12
