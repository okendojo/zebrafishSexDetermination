#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=raconPolishing

#load the modules

module load winnowmap
module load samtools
module load meryl
module add minimap2
module add merfin
module add bcftools
module add racon
module add pbipa 


cd /data/okendojo/zebrafish/data/fish11/polishing/racon

#Run racon polisher
racon -t 24 -q 20 illumina_reads.fastq.gz alignment.sam fish11.fasta > polished_fish11.fasta

