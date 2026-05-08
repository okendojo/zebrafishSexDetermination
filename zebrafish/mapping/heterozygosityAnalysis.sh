#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=heterozygosity

#load the modules
module add minimap2
module add bedtools

cd /data/okendojo/zebrafish/data/fish11/gggenome

# MF2v0.8. maske cent
mat_fa=/data/okendojo/zebrafish/data/fish11/gggenome/ref.fasta
pat_fa=/data/okendojo/zebrafish/data/fish11/gggenome/fish11_merfin_polished.fasta

# split to chr
fastakit -ap mat -o MF2_mat.v0.8.addTag.fasta $mat_fa
fastakit -ap pat -o MF2_pat.v0.8.addTag.fasta $pat_fa

python3 CN1/bin/split_genome_bychrs.py MF2_mat.v0.8.addTag.fasta MF2_pat.v0.8.addTag.fasta

# autosome
for i in `seq 1 22`
do
 chr="chr"$i
 [ -d $chr ] || mkdir $chr
 cd $chr
 echo "#!/bin/bash
date
minimap2 -x asm5  -t 24 --cs mat_${chr}.fa pat_${chr}.fa |sort -k6,6 -k8,8n  > $chr.minimap.sort.paf
paftools.js call $chr.minimap.sort.paf > $chr.var.txt
date" >  $chr.minimap.sh
 #sbatch -c 24 --mem 16g $chr.minimap.sh
 grep -v '^R' $chr.var.txt |cut -f 2-4|sed 's/mat_//g' > $chr.var.bed
 sort -k1,1V -k2n $chr.var.bed |bedtools merge -i - > $chr.var.bed.merge
 cd ..
done

samtools faidx MF2_mat.v0.8.addTag.fasta

cut -f 1,2 MF2_mat.v0.8.addTag.fasta.fai > mat.v0.8.genome.size

bedtools makewindows -g mat.v0.8.genome.size  -w 500000 -s 500000 > mat.genome.500k.bed

### all kinds of variants
cat chr*/*.var.txt |grep -v '^R' >autosome.allvar.txt
cat chr*/*.var.bed.merge > autosome.allvar.bed
bedtools coverage -a ../2.collect/mat.genome.500k.bed -b autosome.allvar.bed > autosome.500k.allvar.cov

# h=0.0004 is the threshold to assign block to homo-(single path with grey color) or heter- (doble paths with other colors)
python3 CN1/bin/makeGFA_from_VARcov.py autosome.500k.allvar.cov 0.0004 > autosome.500k.allvar.0.0004.ab.gfa

# then visualize autosome.500k.allvar.0.0004.ab.gfa using Bandage

