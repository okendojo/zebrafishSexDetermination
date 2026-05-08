#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=RagTag_scaffold

#load the modules

module load ragtag

#move to the right dir
cd /data/okendojo/alono/bacterial_eqa/assembly/Bact_002

ragtag.py scaffold Shigella_sonnei.fasta bac002.fasta  -w -t 24 
