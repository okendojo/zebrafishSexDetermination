#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=600g
#SBATCH --ntasks-per-core=1
#SBATCH --time=25:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=blastn

#Load the module
module load blat

cd /data/okendojo/zebrafish/data/AB/polishing/rDNAcopynumber

#blat ../../asm_verkko2_2/assembly.fasta /data/okendojo/zebrafish/data/fish6/asm/tangleResolution/rDNA.fasta  out.psl

blat -minIdentity=95 /data/okendojo/zebrafish/data/fish6/ontData/concatenated/ont_filt.fasta rdna.fa F6_blat.psl
blat -minIdentity=95 /data/okendojo/zebrafish/data/fish11/ont/concatenate/ont.fasta  rdna.fa F11_blat.psl
blat -minIdentity=95 /data/okendojo/zebrafish/data/AB/batches_ont/ont_ul.fasta rdna.fa AB_blat.psl

#blat Xt_notch4.fasta liver_reads.fasta -t=dna -q=dna Xtmp.psl
