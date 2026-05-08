#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=meryl_db_creation

# Load meryl module
module add meryl

#move the directory containing the data
cd /data/okendojo/zebrafish/data/g3/illumina/meryldbs


# AB only
meryl output AB.only.meryl difference [ difference [ difference AB.meryl TL.meryl ] TU.meryl ] WIK.meryl

# TL only
meryl output TL.only.meryl difference [ difference [ difference TL.meryl AB.meryl ] TU.meryl ] WIK.meryl

# TU only
meryl output TU.only.meryl difference [ difference [ difference TU.meryl WIK.meryl ] AB.meryl ] TL.meryl

# WIK only
meryl output WIK.only.meryl difference [ difference [ difference WIK.meryl TU.meryl ] AB.meryl ] TL.meryl
