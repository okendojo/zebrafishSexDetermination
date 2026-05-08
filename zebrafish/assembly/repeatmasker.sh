#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=repmasking


#Load the modules
module load repeatmasker

cd /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/centromere

RepeatMasker -e ncbi -pa $(expr $SLURM_CPUS_PER_TASK / 2)  -xsmall -gff -dir cent2_mask -species "vertebrates" -html /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/asm_polished/NHGRI_Fish6_cons.fasta
