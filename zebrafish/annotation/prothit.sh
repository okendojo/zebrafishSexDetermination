#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=170:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=protHint

#Load the required modules
module add repeatmasker


cd /data/okendojo/zebrafish/data/fish11/annotation

prothint.py fish11_merfin_polished.fasta uniprotkb_zebrafish_2023_10_05.fasta --workdir protHint
