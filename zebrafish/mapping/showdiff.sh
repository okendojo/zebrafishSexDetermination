#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=nucmer

#load the modules
module load mummer

cd /data/okendojo/zebrafish/data/fish6/asm/variationAnalysis

delta-filter -1 fDanRef_fish11.delta > fDanRef_fish11_filtered.delta

delta-filter -1 GRCz11_fish11.delta > GRCz11_fish11_filtered.delta

show-diff -r fDanRef_fish11_filtered.delta > fDanRef_fish11_filtered_diff.delta

show-diff -r GRCz11_fish11_filtered.delta > GRCz11_fish11_filtered_diff.delta
