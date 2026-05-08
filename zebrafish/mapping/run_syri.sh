#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=syri

#load the modules

module load minimap2
module load syri
module load plotsr

cd /data/okendojo/zebrafish/data/fish6/asm/variationAnalysis


#run syri
syri -c GRCz11_fish11.sam -r chr1_25.fasta -q fish11.fasta --prefix GRCz11_fish11 -F S --nc 24
