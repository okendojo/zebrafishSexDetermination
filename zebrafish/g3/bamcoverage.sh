#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=bamcoverage



#Load the module
module load deeptools
module load samtools

#cd /data/okendojo/zebrafish/data/g3/illumina/bamCoverage

cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/ncb_ref_analysis

#AB="/data/okendojo/zebrafish/data/g3/illumina/AB_F_CB/raw/AB.dedup.bam"
#TL="/data/okendojo/zebrafish/data/g3/illumina/TL_F_CB/raw/TL.dedup.bam"
#TU="/data/okendojo/zebrafish/data/g3/illumina/TU_M_CB/raw/TU.dedup.bam"
#WIK="/data/okendojo/zebrafish/data/g3/illumina/WIK_M_CB/raw/WIK.dedup.bam"

for file in *.bam ; do samtools index ${file} ; done # Index bam files

for bam in *.bam ; do bamCoverage -b $bam --outFileName ${bam}.bw --outFileFormat bigwig --verbose ; done 



#bamCoverage -b ${AB} --outFileName AB.bw --outFileFormat bigwig --verbose 

#bamCoverage -b ${TL} --outFileName TL.bw --outFileFormat bigwig --verbose

#bamCoverage -b ${TU} --outFileName TU.bw --outFileFormat bigwig --verbose

#bamCoverage -b ${WIK} --outFileName WIK.bw --outFileFormat bigwig --verbose





