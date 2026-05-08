#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=vg_mpmap_autoindex

#load modules

module load pangenome


cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/vcf_data/completeData

mkdir -p  vgmapTempDir


#Run the vg autoindex files
vg autoindex --threads 24 --workflow mpmap --workflow rpvg --prefix F0_quant --ref-fasta chr1_25.fasta --vcf f0_phasedDNAvariants.vcf.gz --tx-gff genes.gtf 


