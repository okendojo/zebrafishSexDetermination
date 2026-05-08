#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=chr4

#load module
module add minimap2

#Set data paths
GRCz11=/vf/users/okendojo/zebrafish/data/fish11/mapping/chr4_grcz11.fasta
FISH11=/vf/users/okendojo/zebrafish/data/fish11/mapping/chr4_t2t.fasta

#move to the working directory
cd /data/okendojo/zebrafish/data/fish11/mapping/minimap2

minimap2 -x asm5 -t 24 -o chr4_alignment.paf -L -k 21 $GRCz11 $FISH11 
