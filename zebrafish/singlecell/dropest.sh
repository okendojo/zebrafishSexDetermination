#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=dropest

#Load the module
module use --append ~/modulefiles
module load dropest

cd /data/okendojo/scRNA_project/hui_data/shawnRerun

BAM=/data/okendojo/scRNA_project/hui_data/shawnRerun/24ARC_count_output/24regen_multiome/outs/gex_possorted_bam.bam
CONFIG=/data/okendojo/scRNA_project/hui_data/shawnRerun/dropEst/configs/10x.xml 

dropest  -c ${CONFIG} -o 24hpa -R  -b ${BAM}
 
