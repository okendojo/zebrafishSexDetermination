#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=10:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=synteny

#load the modules

module load minimap2
module load syri
module load plotsr
module load samtools
module load sniffles

cd /data/okendojo/zebrafish/data/assembliesVariants

minimap2 -ax map-hifi -t 30 -k 21 --eqx /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/asm_polished/NHGRI_Fish6_cons.fasta  /data/okendojo/zebrafish/data/fish6/hifi/hifi.fasta -a -o GRCz12tu_ab.sam
samtools sort -O BAM -o GRCz12tu_ab.bam --threads 21  GRCz12tu_ab.sam 
samtools index GRCz12tu_ab.bam

sniffles --input GRCz12tu_ab.bam --vcf GRCz12tu_ab.vcf --reference /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/asm_polished/NHGRI_Fish6_cons.fasta --threads 21 


minimap2 -ax map-hifi -t 30 -k 21 --eqx /data/okendojo/zebrafish/data/fish11/finalPolish/nhgri_polished_fish11.fasta /data/okendojo/zebrafish/data/fish11/hifi/hifi.fasta -a -o GRCz12tu_f11.sam
samtools sort -O BAM -o GRCz12tu_f11.bam --threads 21  GRCz12tu_f11.sam
samtools index GRCz12tu_f11.bam

sniffles --input GRCz12tu_f11.bam --vcf GRCz12tu_f11.vcf --reference /data/okendojo/zebrafish/data/fish11/finalPolish/nhgri_polished_fish11.fasta --threads 21 

minimap2 -ax map-hifi -t 30 -k 21 --eqx /data/okendojo/zebrafish/data/AB/polishing/asm/GRCz12ab.fasta /data/okendojo/zebrafish/data/AB/hifi/*.fastq.gz -a -o GRCz12ab.sam
samtools sort -O BAM -o GRCz12ab.bam --threads 21  GRCz12ab.sam
samtools index GRCz12ab.bam

sniffles --input GRCz12ab.bam --vcf GRCz12ab.vcf --reference /data/okendojo/zebrafish/data/AB/polishing/asm/GRCz12ab.fasta --threads 21 

#run syri
#syri -c GRCz12ab_Dr1AB.bam  --prefix GRCz12ab_Dr1AB.bam -r GRCz12ab.fasta  -q Dr1AB.fasta  -k -F B --nc 24

#plotsr --sr GRCz12ab_Dr1AB.bamsyri.out --genome genomes.txt -o GRCz12ab_Dr1AB1.pdf -R -H 8 -W 5 -b pdf --chrord chrs2.txt

