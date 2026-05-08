#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=120g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=compare_asm

#load the modules

module load winnowmap
module load samtools
module load python

#move to the right dir
cd /data/okendojo/zebrafish/data/fish11/asm_mapping

python compareAssemblies/compareAssemblies.py -r verkko_assembly.fasta -q hifiasm_assembly.fasta -o compare_asm -n 10
