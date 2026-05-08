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

cd  /data/okendojo/paradisfishProject/centromere

RepeatMasker -e ncbi -pa $(expr $SLURM_CPUS_PER_TASK / 2)  -dir centromereseqs_mask -species "vertebrates" cent.fasta
