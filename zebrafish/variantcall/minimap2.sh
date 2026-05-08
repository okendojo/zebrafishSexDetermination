#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=64g
#SBATCH --ntasks-per-core=1
#SBATCH --time=96:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov

#load the module
module load bwa
module load minimap2
module load samtools
module load dragmap/1.2.1
module load mummer

cd /data/okendojo/paradisfishProject/genoComparison/

minimap2 -ax asm5 -t 16 -k 21 --eqx betta2.fasta macope.fasta -o macOpe2.sam

samtools view -S -b macOpe2.sam -o macOpe2.bam

samtools sort macOpe2.bam -o macOpe2.sorted.bam

samtools index macOpe2.sorted.bam

#mummer -maxmatch -l 20 -b -n -k 3 -threads 16 Betta_splendens.fa mc.fasta




