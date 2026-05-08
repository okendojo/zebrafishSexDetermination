#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=mkref

#Load the module
module load cellranger
module load cellranger-arc

cd /data/okendojo/scRNA_project/hui_data/rawData/fastqs_ATAC

cellranger-arc mkref --config=/home/okendojo/scripts/zebrafish/singlecell/grcz11.config --nthreads=24 --memgb=200

