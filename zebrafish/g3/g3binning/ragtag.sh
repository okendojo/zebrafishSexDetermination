#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=ragtag_scaffold

#load modules
module load ragtag
module load seqkit

cd /data/okendojo/zebrafish/data/g3/assembly/tu/tu_asm

ragtag.py scaffold  /data/okendojo/zebrafish/data/fish6/asm/variationAnalysis/danRer4.fasta  assembly.fasta -o ragtag -t 24 -w -C  -f 1000 --remove-small

cd /data/okendojo/zebrafish/data/g3/assembly/tl/tl_asm
ragtag.py scaffold  /data/okendojo/zebrafish/data/fish6/asm/variationAnalysis/danRer4.fasta  assembly.fasta -o ragtag -t 24 -w -C -f 1000 --remove-small

cd /data/okendojo/zebrafish/data/g3/assembly/wik/wik_asm

ragtag.py scaffold  /data/okendojo/zebrafish/data/fish6/asm/variationAnalysis/danRer4.fasta  assembly.fasta -o ragtag -t 24 -w -C -f 1000 --remove-small

cd /data/okendojo/zebrafish/data/g3/assembly/ab/ab_asm

ragtag.py scaffold  /data/okendojo/zebrafish/data/fish6/asm/variationAnalysis/danRer4.fasta  assembly.fasta -o ragtag -t 24 -w -C  -f 1000 --remove-small 



