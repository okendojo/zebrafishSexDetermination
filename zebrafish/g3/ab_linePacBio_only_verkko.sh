#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=abLine_PacBio_only


#Load the modules
module load verkko/2.1


#Run the reads extraction first

cd /data/okendojo/zebrafish/data/g3/assembly

#Run verkko
verkko -d  abLine_PacBio_only  --hifi /data/okendojo/zebrafish/data/g3/ab_inbreadPacBio/*.fastq.gz   --no-nano  --slurm --lay-run 1 64 24
