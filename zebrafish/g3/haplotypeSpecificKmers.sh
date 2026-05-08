#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=90:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=haplotypeKmers

#load modules
module add meryl

cd /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls


# AB only
meryl output AB.only.meryl difference [ difference [ difference AB.k21.meryl TL.k21.meryl ] TU.k21.meryl ] WIK.k21.meryl

# TL only
meryl output TL.only.meryl difference [ difference [ difference TL.k21.meryl AB.k21.meryl ] TU.k21.meryl ] WIK.k21.meryl

# TU only
meryl output TU.only.meryl difference [ difference [ difference TU.k21.meryl WIK.k21.meryl ] AB.k21.meryl ] TL.k21.meryl

# WIK only
meryl output WIK.only.meryl difference [ difference [ difference WIK.k21.meryl TU.k21.meryl ] AB.k21.meryl ] TL.k21.meryl
