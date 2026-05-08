#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=64g
#SBATCH --ntasks-per-core=1
#SBATCH --time=96:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov

#load module
module load smoove
module load manta

#cd /data/okendojo/paradisfishProject/snv/dragmapping

REF=/data/okendojo/datashare/macOpeProject/macOpe2Assembly.fasta

cd /data/okendojo/paradisfishProject/snv/dragmapping/mantra


