#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=170:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=BLAT

#Load the required modules
module add blat


cd /data/okendojo/zebrafish/data/fish11/annotation

blat fish11_merfin_polished.fasta centromereproteins.fasta  -t=dna  -q=dna fish11_centpaPc_gene.psl -out=psl
