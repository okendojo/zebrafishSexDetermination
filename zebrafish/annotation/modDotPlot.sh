#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=ModDotPlot


#move to the right directory

DIR="/data/okendojo/zebrafish/data/fish11/gggenome/original"

cd ${DIR}


moddotplot -i fish11_merfin_polished.fasta -id 80 -o moddotPlots
