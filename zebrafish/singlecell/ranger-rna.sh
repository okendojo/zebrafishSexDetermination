#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=ranger-rna

#Load the module
module load cellranger

cd /data/okendojo/scRNA_project/hui_data/rawData/fastqs_GEX

#Run 24hpa
cellranger count --id=run_count_24hpa \
   --fastqs=24hpa \
   --sample=24hpa_S1 \
   --transcriptome=/vf/users/okendojo/scRNA_project/hui_data/rawData/fastqs_GEX/GRCz11 --chemistry ARC-v1	

#Run 36hpa
cellranger count --id=run_count_36hpa \
   --fastqs=36hpa \
   --sample=36hpa_S3 \
   --transcriptome=/vf/users/okendojo/scRNA_project/hui_data/rawData/fastqs_GEX/GRCz11 --chemistry ARC-v1 

#Run 48 hpa
cellranger count --id=run_count_48hpa \
   --fastqs=48hpa \
   --sample=48hpa_S4 \
   --transcriptome=/vf/users/okendojo/scRNA_project/hui_data/rawData/fastqs_GEX/GRCz11 --chemistry ARC-v1

#Run UI_1
cellranger count --id=run_count_UI_1 \
   --fastqs=UI_1 \
   --sample=UI_S2 \
   --transcriptome=/vf/users/okendojo/scRNA_project/hui_data/rawData/fastqs_GEX/GRCz11 --chemistry ARC-v1








