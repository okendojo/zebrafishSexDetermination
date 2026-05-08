#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=520g
#SBATCH --ntasks-per-core=1
#SBATCH --time=40:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=fetchngs

#load modules
module add nextflow
module add singularity


#cd /data/okendojo/zebrafish/data/g3/sex_project/pintoanalysis/fetchngs #pinto data

cd /data/okendojo/zebrafish/data/g3/sex_project/asm/t2tMap/bam_fixed_rg

#nextflow run nf-core/fetchngs -profile biowulf --input accession.csv  -resume  --outdir fastq  -c  custom.config   --download_method sratools   --email 'javanokendo@gmail.com'

nextflow run nf-core/fetchngs --input SRR_Acc_List.tsv  -resume  -profile biowulf --outdir na_rna_fastq

