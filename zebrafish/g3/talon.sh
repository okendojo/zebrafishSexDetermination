#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=210g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=isoseq

#load modules
module load talon

cd /data/okendojo/zebrafish/data/annotation/fish6/annot/

talon --f completeconfig.csv --db talon.db --o grcz12tuannot --verbosity 6 --threads 24 --build /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/asm_polished/NHGRI_Fish6_cons.fasta
