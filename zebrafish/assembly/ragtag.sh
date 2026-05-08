#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=230g
#SBATCH --ntasks-per-core=1
#SBATCH --time=10:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=ragtag_scaffold

#load modules
module load ragtag
module load seqkit

cd  /data/okendojo/zebrafish/data/g3/assembly/asm_270625

ragtag.py scaffold  /data/ShawnBurgessLab/javan_duncan_multiome/grcz12tuannotations/GRCz12tu.fasta assembly.fasta -o g3_scaffold -q 20 -C -r -w -g 10 -t 24 -w
#ragtag.py scaffold  /data/ShawnBurgessLab/javan_duncan_multiome/grcz12tuannotations/GRCz12tu.fasta ZW.hap2.fasta -o ZW_hap2 -q 20 -C -r -w -g 10 -t 24 -w
#ragtag.py scaffold  /data/ShawnBurgessLab/javan_duncan_multiome/grcz12tuannotations/GRCz12tu.fasta /vf/users/okendojo/zebrafish/data/g3/sex_project/scinkd_analysis/8WW.hap1.fasta -o 8WW_hap1 -q 20 -C -r -w -g 10 -t 24 -w
#ragtag.py scaffold  /data/ShawnBurgessLab/javan_duncan_multiome/grcz12tuannotations/GRCz12tu.fasta /vf/users/okendojo/zebrafish/data/g3/sex_project/scinkd_analysis/8WW.hap2.fasta -o 8WW_hap2 -q 20 -C -r -w -g 10 -t 24 -w
#ragtag.py scaffold  /data/ShawnBurgessLab/javan_duncan_multiome/grcz12tuannotations/GRCz12tu.fasta /vf/users/okendojo/zebrafish/data/g3/sex_project/scinkd_analysis/CBZW.hap1.fasta -o CBZW_hap1 -q 20 -C -r -w -g 10 -t 24 -w
#ragtag.py scaffold  /data/ShawnBurgessLab/javan_duncan_multiome/grcz12tuannotations/GRCz12tu.fasta /vf/users/okendojo/zebrafish/data/g3/sex_project/scinkd_analysis/CBZW.hap2.fasta -o CBZW_hap2 -q 20 -C -r -w -g 10 -t 24 -w


#ragtag.py scaffold  /data/ShawnBurgessLab/javan_duncan_multiome/grcz12tuannotations/GRCz12tu.fasta W.hap1.fasta.gz -o hap1 -q 20 -C -r -w -g 10 -t 24 -w
#ragtag.py scaffold  /data/ShawnBurgessLab/javan_duncan_multiome/grcz12tuannotations/GRCz12tu.fasta W.hap2.fasta.gz -o hap2 -q 20 -C -r -w -g 10 -t 24 -w

