#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=maker_001_annotZF

#Load the required modules
#module add repeatmasker
module add maker/3.01.03
module load augustus
module load busco
module load blast
module load genometools
module load seqkit
module load bioawk

#======================
###paradisefish genome annotation##
#=====================
#1. Run maker annotation

cd /data/okendojo/zebrafish/data/fish11/annotation/annot_01_maker

maker -base maker_01 maker_opts.ctl maker_bopts.ctl maker_exe.ctl 
