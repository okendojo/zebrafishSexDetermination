#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=centrominer


#Load modules
module load python/3.9 minimap2 trf cd-hit mummer/4.0.0beta2  R gnuplot/5.4.3 blast


cd  /data/okendojo/paradisfishProject/centromere

#1 Fish 6
#te="/data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/centromere/NHGRI_Fish6_cons.fasta.mod.EDTA.TEanno.gff3"
#gene="/data/okendojo/zebrafish/data/fish6/asm/t2t/liftoff/NHGRI_Fish6.gtf"

#quartet_centrominer.py -i /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/asm_polished/NHGRI_Fish6_cons.fasta --TE ${te} --gene ${gene} -t 24 -p NHGRI_Fish6_CEN  -s 0.99 -l 1000000

#2 Fish11
quartet_centrominer.py -i t2t_macope2.fasta -p MacOpe_CEN -t 24 -s 0.90 
