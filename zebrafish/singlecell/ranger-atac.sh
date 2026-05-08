#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=ranger-atac

#Load the module
module load cellranger-atac

#move to the working directory
cd /data/okendojo/scRNA_project/hui_data/rawData/fastqs_ATAC

cellranger-atac count --id=run_count_24hpa \
   --fastqs=/data/okendojo/scRNA_project/hui_data/rawData/fastqs_ATAC/24hpa \
   --sample=24hpa \
   --reference=/data/okendojo/scRNA_project/hui_data/rawData/fastqs_ATAC/GRCz11 \
   --jobmode=slurm --noexit 


cellranger-atac count --id=run_count_36hpa \
   --fastqs=/data/okendojo/scRNA_project/hui_data/rawData/fastqs_ATAC/36hpa \
   --sample=36hpa \
   --reference=/data/okendojo/scRNA_project/hui_data/rawData/fastqs_ATAC/GRCz11 \
   --jobmode=slurm --noexit


cellranger-atac count --id=run_count_48hpa \
   --fastqs=/data/okendojo/scRNA_project/hui_data/rawData/fastqs_ATAC/48hpa \
   --sample=48hpa \
   --reference=/data/okendojo/scRNA_project/hui_data/rawData/fastqs_ATAC/GRCz11 \
   --jobmode=slurm --noexit

cellranger-atac count --id=run_count_UI_1 \
   --fastqs=/data/okendojo/scRNA_project/hui_data/rawData/fastqs_ATAC/UI_1 \
   --sample=UI \
   --reference=/data/okendojo/scRNA_project/hui_data/rawData/fastqs_ATAC/GRCz11 \
   --jobmode=slurm --noexit
