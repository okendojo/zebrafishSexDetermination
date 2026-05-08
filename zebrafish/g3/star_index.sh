#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=100g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=star_index

#load modules
module add nextflow
module add singularity
module add STAR

cd  /data/okendojo/zebrafish/data/g3/F2_variants

ref="/data/okendojo/zebrafish/data/fish6/asm/t2t/NHGRI_Fish6.fasta"
gtf="/data/okendojo/zebrafish/data/fish6/liftoff/fish6_asm/fish6_annotation.gtf"
#gtf="/data/okendojo/zebrafish/data/g3/F2_variants/fish6_fixedAgatfixed.gtf"

STAR --runMode genomeGenerate --sjdbOverhang 149 --genomeDir  grcz12_star_index --genomeFastaFiles $ref --sjdbGTFfile  $gtf


#refwik="wiktl.fasta"
#gtfwik="wiktl.gtf"

#STAR --runMode genomeGenerate --sjdbOverhang 149 --genomeDir  wiktu_star_index --genomeFastaFiles $refwik --sjdbGTFfile  $gtfwik
