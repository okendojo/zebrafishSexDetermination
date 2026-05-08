#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=polishing

#load modules
module load pbipa
module load samtools bcftools winnowmap merfin meryl racon

cd /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/variantCall

bash automated-polishing.sh 24 2 /data/okendojo/zebrafish/data/fish6/asm/t2t/NHGRI_Fish6.fasta  /data/okendojo/zebrafish/data/fish6/hifi/m54313U_220817_024630.hifi_reads.fastq.gz /data/okendojo/zebrafish/data/fish6/asm/polishingasm/hifi_meryl/PacBio.k21.meryl hifi_polishing
