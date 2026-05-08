#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=polishing

#load the modules

module load winnowmap
module load samtools
module load meryl
module add minimap2
module add merfin
module add bcftools
module add racon
module add pbipa 

cd /data/okendojo/zebrafish/data/fish11/polishing

sh /data/okendojo/zebrafish/data/fish11/mapping/T2T-Polish/automated_polishing/automated-polishing.sh 24 2 /data/Zebrafish_T2T/fish11/resolved_assembly/fish11.fasta /data/okendojo/zebrafish/data/fish11/hifi/hifi.fasta /data/okendojo/zebrafish/data/fish11/polishing/repetitive_k21.txt HiFi_polish
