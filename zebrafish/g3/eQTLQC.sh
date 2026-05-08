#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=200g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=eQTLQC

#load modules
module add GATK
module load R  plink/1.9.0-beta4.4 
module add module load bowtie2
module load rsem

cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/qtlTMP_analysis/eQTLQC/Sample/fastq

./run_demo.sh
