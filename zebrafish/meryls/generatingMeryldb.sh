#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=80g
#SBATCH --ntasks-per-core=1
#SBATCH --time=96:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#move to the working directory

#Load the modules
module load merqury/1.3

cd /data/okendojo/zebrafish/data/g3/illumina/meryldbs

meryl count k=21 /data/okendojo/zebrafish/data/g3/illumina/AB_F_CB/qc/*.gz output AB.meryl
meryl count k=21 /data/okendojo/zebrafish/data/g3/illumina/TL_F_CB/qc/*.gz output TL.meryl
meryl count k=21 /data/okendojo/zebrafish/data/g3/illumina/TU_M_CB/qc/*.gz output TU.meryl
meryl count k=21 /data/okendojo/zebrafish/data/g3/illumina/WIK_M_CB/qc/*.gz output WIK.meryl
meryl count k=21 /data/okendojo/zebrafish/data/g3/ontData/individual/*.fastq.gz output child.meryl

