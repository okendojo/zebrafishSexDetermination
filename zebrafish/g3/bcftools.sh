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
module load bcftools

cd  /data/okendojo/zebrafish/data/segmental_duplication/synteny/SV

#SET PATHS
Ref_FASTA="/data/okendojo/zebrafish/data/g3/illumina/WIK_M_CB/raw/Danio_rerio.GRCz11.dna.primary_assembly.fa"
WIK="/data/okendojo/zebrafish/data/g3/illumina/WIK_M_CB/raw/bamFile.txt"
AB="/data/okendojo/zebrafish/data/g3/illumina/AB_F_CB/raw/filebam.txt"
TL="/data/okendojo/zebrafish/data/g3/illumina/TL_F_CB/raw/file.txt"
TU="/data/okendojo/zebrafish/data/g3/illumina/TU_M_CB/raw/bamfile.txt "


# 1, call variants
bcftools mpileup -Ou --redo-BAQ --min-BQ 20   --threads 24   -f /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/asm_polished/NHGRI_Fish6_cons.fasta -b fileList.txt | bcftools call -mv -Oz --threads 24 -o f6f11_snps.vcf.gz

#bcftools mpileup -Ou --redo-BAQ --min-BQ 20   --threads 24 --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f ${Ref_FASTA} -b ${AB} | bcftools call -mv -Oz --threads 24 -o AB.vcf.gz

#bcftools mpileup -Ou --redo-BAQ --min-BQ 20   --threads 24 --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f ${Ref_FASTA} -b ${TL} | bcftools call -mv -Oz --threads 24 -o TL.vcf.gz

#bcftools mpileup -Ou --redo-BAQ --min-BQ 20   --threads 24 --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f ${Ref_FASTA} -b ${TU} | bcftools call -mv -Oz --threads 24 -o TU.vcf.gz
