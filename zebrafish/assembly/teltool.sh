#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=teltool2


module load samtools

cd /data/okendojo/zebrafish/data/fish6/centromere/teltools

#perl telsize/extractTelReads/extract_telomeric_reads.pl /data/okendojo/zebrafish/data/fish6/centromere/teltools/fulTel.fasta fish6_telomeres

#python telsize/estimateTelLength/estimateTelSize.py  --motif CCTAAC --noseq fish6_telomeres.telomeric.fasta.gz > f6_telsize.txt 

samtools view -@24 -O BAM -o telomere_alignment.bam -bS telomere_alignment.sam

python telfuse/extract_sites_sample.py telomere_alignment.bam fish6.fasta fish6_telfuse 
