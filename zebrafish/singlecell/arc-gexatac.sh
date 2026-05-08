#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=mkref-rna

#Load the module
module load cellranger
module load cellranger-arc

cd /data/okendojo/scRNA_project/hui_data/rawData
# 24 hpa
#cellranger-arc count --id=24hpa \
#                       --reference=/data/okendojo/scRNA_project/hui_data/rawData/fastqs_ATAC/GRCz11 \
 #                      --libraries=/data/okendojo/scRNA_project/hui_data/rawData/24hpa.csv \
  #                     --localcores=16 \
   #                    --localmem=64

# 36 hpa
#cellranger-arc count --id=36hpa \
 #                      --reference=/data/okendojo/scRNA_project/hui_data/rawData/fastqs_ATAC/GRCz11 \
  #                     --libraries=/data/okendojo/scRNA_project/hui_data/rawData/36hpa.csv \
   #                    --localcores=16 \
    #                   --localmem=64

# 48 hpa
#cellranger-arc count --id=48hpa \
 #                      --reference=/data/okendojo/scRNA_project/hui_data/rawData/fastqs_ATAC/GRCz11 \
  #                     --libraries=/data/okendojo/scRNA_project/hui_data/rawData/48hpa.csv \
   #                    --localcores=16 \
    #                   --localmem=64

# UI
cellranger-arc count --id=UI \
                       --reference=/data/okendojo/scRNA_project/hui_data/rawData/fastqs_ATAC/GRCz11 \
                       --libraries=/data/okendojo/scRNA_project/hui_data/rawData/UI.csv \
                       --localcores=16 \
                       --localmem=64








