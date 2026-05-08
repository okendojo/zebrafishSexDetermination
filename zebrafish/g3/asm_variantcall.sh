#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=120g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=bcftools_variantCall

#load module
module load minimap2
module load samtools
module load k8

cd /data/okendojo/zebrafish/data/segmental_duplication/synteny/SV

minimap2 -x asm5 -c --cs=long /data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.dna.primary_assembly.fa /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/asm_polished/NHGRI_Fish6_cons.fasta > GRCz11_GRcz12tu.paf
sort -k6,6 -k8,8n GRCz11_GRcz12tu.paf | paftools.js call - > GRCz11_GRCz12tu_variants.vcf


#samtools view -b alignment.sam | samtools sort -o alignment.bam
#samtools index alignment.bam

