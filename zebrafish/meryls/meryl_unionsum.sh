#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=96:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov

cd /data/okendojo/zebrafish/data/g3/illumina/meryldbs

meryl union-sum AB.meryl/ TU.meryl/ output tu_ab.meryl

meryl union-sum WIK.meryl TL.meryl output wik_tl.meryl
