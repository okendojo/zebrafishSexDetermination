#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=96:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov

cd /data/okendojo/zebrafish/data/g3/illumina/meryldbs

# AB only
meryl  difference [ difference [ difference AB.meryl TL.meryl ] TU.meryl ] WIK.meryl output AB.only.meryl

# TL only
meryl difference [ difference [ difference TL.meryl AB.meryl ] TU.meryl ] WIK.meryl output TL.only.meryl

# TU only
meryl difference [ difference [ difference TU.meryl WIK.meryl ] AB.meryl ] TL.meryl output TU.only.meryl

# WIK only
meryl difference [ difference [ difference WIK.meryl TU.meryl ] AB.meryl ] TL.meryl output WIK.only.meryl
