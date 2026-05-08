#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=telowinnowmap

#load the modules

module load winnowmap
module load samtools

#move to the right dir
cd /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/pur


#count the kmers

##meryl count k=21 output meryldb fish11_assembly.fasta

##meryl print greater-than distinct=0.9998 meryldb > repetitive_k21.txt

REF="NHGRI_Fish6_cons.fasta"
ASM="GRCz11.fasta"

winnowmap -ax asm20 -t 20 -H --MD -k 21 $REF $ASM > NHGRI_fish6_GRCZ11_mapping.sam 

#winnowmap -x asm5 ${REF} ${ASM}  -k 21 -a  --eqx -t 32 -o NHGRI_fish6_GRCZ11_mapping.sam
