#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=800g
#SBATCH --gres=lscratch:500
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=F11_F6_winnowmap

#load the modules

module load winnowmap
module load samtools

export TMPDIR=/lscratch/$SLURM_JOB_ID

#move to the right dir
cd /data/okendojo/zebrafish/data/fish11/asm_mapping/fish11_fish6_hifiMapping/fish11

#count the kmers
meryl count k=21 output meryldb assembly.fasta

meryl print greater-than distinct=0.9998 meryldb > repetitive_k21.txt

echo "now mapping Fish11 data........."

winnowmap -W repetitive_k21.txt -t 32 -ax map-pb assembly.fasta m64467e_230509_164640.hifi_reads.fastq.gz  m64467e_230511_034148.hifi_reads.fastq.gz > PB_verkko.sam

echo "Now running samtools on ONT data"

samtools view -S -b PB_verkko.sam -o PB_verkko.bam
samtools sort PB_verkko.bam -o PB_verkko-sorted.bam
samtools index PB_verkko-sorted.bam

cd /data/okendojo/zebrafish/data/fish11/asm_mapping/fish11_fish6_hifiMapping/fish6

echo "mapping pacbio on verkko"

#count the kmers
meryl count k=21 output meryldb assembly.fasta

meryl print greater-than distinct=0.9998 meryldb > repetitive_k21.txt

winnowmap -W repetitive_k21.txt -t 32 -ax map-pb assembly.fasta m54313U_220817_024630.hifi_reads.fastq.gz > PB_verkko.sam

echo "running samtools on verkko assembly"
samtools view -S -b PB_verkko.sam -o PB_verkko.bam
samtools sort PB_verkko.bam -o PB_verkko-sorted.bam
samtools index PB_verkko-sorted.bam














