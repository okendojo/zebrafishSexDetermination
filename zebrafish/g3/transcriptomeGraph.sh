#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=90:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=transcriptomeGraph

#Load the modules
module add pangenome

#move to working dir

cd /data/okendojo/zebrafish/data/g3/eQTL

vg rna -t 24 -p --transcripts /data/okendojo/zebrafish/refGenome/GRCz11_genomic.gtf -u -b transcriptome_GBTW -g -f transcriptome  -i transcriptome_tsv -r -a 
