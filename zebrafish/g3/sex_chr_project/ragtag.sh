#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=230g
#SBATCH --ntasks-per-core=1
#SBATCH --time=05:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=ragtag_scaffold

#load modules
module load ragtag
module load seqkit minimap2


#process contigs

cd /data/okendojo/zebrafish/data/g3/sex_project/asm

#mkdir -p haplotypes

#for gfa in *.hap1.p_ctg.gfa *.hap2.p_ctg.gfa; do
  # skip loop if no files match
#  [ -e "$gfa" ] || continue

#  base=$(basename "$gfa" .gfa)
#  awk '/^S/{print ">"$2; print $3}' "$gfa" > "haplotypes/${base}.fa"
#done


cd /data/okendojo/zebrafish/data/g3/sex_project/asm/haplotypes

for fa in *.fa; do
  base=$(basename "$fa" .fa)

  ragtag.py scaffold \
     /data/okendojo/zebrafish/data/g3/sex_project/t2t_pangenome/asm/GRCz12tu.fasta \
    "$fa" -w  --remove-small -f 1000 -t 32 \
    -o scaffold/"$base"
done


#Rename ragtags
cd /data/okendojo/zebrafish/data/g3/sex_project/asm/haplotypes/scaffold

for f in *.fa; do
  sed 's/_RagTag//g' "$f" > "${f%.fa}.clean.fa"
done

