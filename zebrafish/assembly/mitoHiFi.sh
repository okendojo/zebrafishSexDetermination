#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=mitoHiFi


cd /data/okendojo/zebrafish/data/mtDNA/f11_mtasm

module load singularity

source /usr/local/current/singularity/app_conf/sing_binds

singularity exec --cleanenv /data/okendojo/zebrafish/data/fish6/asm/polishingasm/finalissues/mito/singularity/mitohifi_master.sif mitohifi.py -r /data/okendojo/zebrafish/data/fish11/hifi/hifi.fasta -f sequence.fasta -g sequence.gb -t 24 -o 2

cd /data/okendojo/zebrafish/data/mtDNA/f6_mtasm

singularity exec --cleanenv /data/okendojo/zebrafish/data/fish6/asm/polishingasm/finalissues/mito/singularity/mitohifi_master.sif mitohifi.py -r /data/okendojo/zebrafish/data/fish6/asm/polishingasm/finalissues/mito/hifi.fasta -f sequence.fasta -g sequence.gb -t 24 -o 2 
