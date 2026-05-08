#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=90:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=binning

#load modules
module add python/3.9 meryl R merqury seqkit

set -euo pipefail

cd /data/okendojo/zebrafish/data/g3 

############################################
# CONFIGURATION
############################################
THREADS=24
K=21
THRESHOLD=0.90

# Input directories
###ILLUMINA_DIR=/data/okendojo/zebrafish/data/g3/illumina
###PACBIO_DIR=/data/okendojo/zebrafish/data/g3/pacBio
###ONT=/data/okendojo/zebrafish/data/g3/ontData/ont_uniq.fastq.gz

OUT=output
###mkdir -p $OUT/{meryl,unique_kmers,PacBioGP,ONTGP,bins/{AB,TU,TL,WIK,ambiguous},plots,tmp}

SAMPLES=(AB TU TL WIK)

############################################
# STEP 1: Locate Illumina reads
############################################
###echo "== Locating Illumina reads =="

#for S in "${SAMPLES[@]}"; do
#    R1=$(ls $ILLUMINA_DIR/$S/*_R1_*.fastq.gz)
#    R2=$(ls $ILLUMINA_DIR/$S/*_R2_*.fastq.gz)
#    echo "$S -> $R1 $R2"
#    eval ${S}_R1=$R1
#    eval ${S}_R2=$R2
#done

############################################
# STEP 2: meryl k-mer counting
############################################
###echo "== Counting k-mers =="

#for S in "${SAMPLES[@]}"; do
#    meryl count k=$K threads=$THREADS \
#        output $OUT/meryl/${S}.meryl \
#        $(eval echo \${${S}_R1}) $(eval echo \${${S}_R2})
#done

############################################
# STEP 3: derive unique k-mers
############################################
###echo "== Generating unique k-mers =="

#for S in "${SAMPLES[@]}"; do
#    OTHERS=()
#    for O in "${SAMPLES[@]}"; do
#        [[ "$O" != "$S" ]] && OTHERS+=("$OUT/meryl/$O.meryl")
#    done
#
#    meryl union-sum output $OUT/meryl/OTHERS_${S}.meryl "${OTHERS[@]}"
#
#    meryl difference output $OUT/unique_kmers/${S}.unique.meryl \
#        $OUT/meryl/${S}.meryl $OUT/meryl/OTHERS_${S}.meryl
#done

############################################
# STEP 4: Merge long reads
############################################
###echo "== Preparing long reads =="

#zcat $PACBIO_DIR/*.fastq.gz > $OUT/tmp/pacbio_all.fastq
#zcat $ONT > $OUT/tmp/ont_all.fastq
#seqkit fq2fa --threads 24 $OUT/tmp/pacbio_all.fastq > $OUT/tmp/pacbio_all.fasta
#seqkit fq2fa --threads 24 $OUT/tmp/ont_all.fastq > $OUT/tmp/ont_all.fasta


############################################
# STEP 5: Identify GP-specific k-mers in long reads
############################################
echo "== Extracting GP-specific kmers from long reads =="

for S in "${SAMPLES[@]}"; do
    meryl-lookup -sequence $OUT/tmp/pacbio_all.fasta \
        -mers $OUT/unique_kmers/${S}.unique.meryl \
        -existence -output $OUT/PacBioGP/${S}.hits.txt

    meryl-lookup -sequence $OUT/tmp/ont_all.fasta \
        -mers $OUT/unique_kmers/${S}.unique.meryl \
        -existence -output $OUT/ONTGP/${S}.hits.txt
done

############################################
# STEP 6: Read assignment (Python)
############################################
echo "== Assigning reads =="

python3 <<EOF
from collections import defaultdict

samples = ["AB","TU","TL","WIK"]
threshold = $THRESHOLD

def assign_reads(prefix):
    counts = defaultdict(lambda: {s:0 for s in samples})

    for s in samples:
        with open(f"$OUT/{prefix}/{s}.hits.txt") as f:
            for line in f:
                rid = line.split()[0]
                counts[rid][s] += 1

    assignments = {}
    for rid, d in counts.items():
        total = sum(d.values())
        if total == 0:
            assignments[rid] = "ambiguous"
            continue
        best = max(d, key=d.get)
        if d[best]/total >= threshold:
            assignments[rid] = best
        else:
            assignments[rid] = "ambiguous"
    return assignments

pb = assign_reads("PacBioGP")
ont = assign_reads("ONTGP")

# Save
with open("$OUT/pacbio_assignments.tsv","w") as f:
    for k,v in pb.items():
        f.write(f"{k}\t{v}\n")

with open("$OUT/ont_assignments.tsv","w") as f:
    for k,v in ont.items():
        f.write(f"{k}\t{v}\n")
EOF

############################################
# STEP 7: Split reads into bins
############################################
echo "== Splitting reads =="

split_fastq () {
    INPUT=$1
    ASSIGN=$2
    PREFIX=$3

    python3 <<EOF
assign = dict(line.strip().split() for line in open("$ASSIGN"))

import os
bins=["AB","TU","TL","WIK","ambiguous"]
files={b:open("$OUT/bins/"+b+"/${PREFIX}.fastq","w") for b in bins}

with open("$INPUT") as f:
    while True:
        h=f.readline().strip()
        if not h: break
        seq=f.readline(); plus=f.readline(); qual=f.readline()
        rid=h.split()[0][1:]
        b=assign.get(rid,"ambiguous")
        files[b].write(h+"\n"+seq+plus+qual)

for f in files.values(): f.close()
EOF
}

split_fastq $OUT/tmp/pacbio_all.fastq $OUT/pacbio_assignments.tsv PacBioReads
split_fastq $OUT/tmp/ont_all.fastq $OUT/ont_assignments.tsv ONTReads

############################################
# STEP 8: Summary stats
############################################
echo "== Generating summary =="

for TYPE in PacBioReads ONTReads; do
    for b in AB TU TL WIK ambiguous; do
        COUNT=$(grep -c "^@" $OUT/bins/$b/${TYPE}.fastq 2>/dev/null || echo 0)
        echo -e "$TYPE\t$b\t$COUNT"
    done
done > $OUT/summary.tsv

############################################
# STEP 9: Merqury QC
############################################
echo "== Running Merqury =="

meryl union-sum output $OUT/meryl/all.meryl $OUT/meryl/*.meryl

# Replace genome.fasta if available
# merqury.sh $OUT/meryl/all.meryl genome.fasta $OUT/merqury

############################################
# STEP 10: R plots
############################################
echo "== Plotting =="

cat << 'RSCRIPT' > $OUT/plot.R
library(ggplot2)
library(readr)
library(dplyr)

df <- read_tsv("summary.tsv", col_names=c("Type","Bin","Reads"))

p1 <- ggplot(df, aes(x=Bin,y=Reads,fill=Bin)) +
  geom_col() + facet_wrap(~Type) +
  theme_bw(base_size=14) +
  theme(legend.position="none")

ggsave("plots/read_distribution.pdf",p1,width=8,height=5)

# Heatmap-style
library(reshape2)
mat <- dcast(df, Type ~ Bin, value.var="Reads")
m <- as.matrix(mat[,-1])
rownames(m) <- mat$Type

library(pheatmap)
pheatmap(m, filename="plots/heatmap.pdf")
RSCRIPT

(cd $OUT && Rscript plot.R)

############################################
echo "Pipeline complete!"
