#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=96:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov

cd /data/okendojo/zebrafish/data/g3/illumina/meryldb

meryl count k=21 /data/okendojo/zebrafish/data/g3/ontData/combined/zont.fastq.gz output child.meryl

$MERQURY/trio/hapmers.sh TU_AB.meryl/ WIK_TL.meryl/ child.meryl


