#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=90:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --gres=lscratch:500
#SBATCH --job-name=contigBinning

#load modules
module add python/3.9 meryl R merqury seqkit

cd /data/okendojo/zebrafish/data/g3/output/unique_kmers

#!/usr/bin/env bash
set -euo pipefail

############################################
# CONFIG
############################################
THREADS=24
THRESHOLD=0.95      # realistic threshold
MIN_KMERS=100        # minimum evidence
ASSEMBLY=/data/okendojo/zebrafish/data/g3/assembly/asm_jan2026/assembly.fasta

SAMPLES=(AB TU TL WIK)

OUT=contig_binning_bed
mkdir -p $OUT/{lookup,results,fasta}

############################################
# STEP 1: meryl lookup with -bed
############################################
echo "== Step 1: meryl lookup (BED mode) =="

for S in "${SAMPLES[@]}"; do
    meryl-lookup \
        -sequence $ASSEMBLY \
        -mers ${S}.unique.meryl \
        -bed \
        -output $OUT/lookup/${S}.bed
done

############################################
# STEP 2: Count k-mer hits per contig
############################################
echo "== Step 2: counting k-mer hits =="

python3 <<EOF
from collections import defaultdict

samples = ["AB","TU","TL","WIK"]
counts = defaultdict(lambda: {s:0 for s in samples})

for s in samples:
    with open(f"$OUT/lookup/{s}.bed") as f:
        for line in f:
            contig = line.split()[0]
            counts[contig][s] += 1

# write table
with open("$OUT/results/contig_counts.tsv","w") as out:
    out.write("contig\tAB\tTU\tTL\tWIK\ttotal\n")
    for c,d in counts.items():
        ab,tu,tl,wik = d["AB"],d["TU"],d["TL"],d["WIK"]
        total = ab+tu+tl+wik
        out.write(f"{c}\t{ab}\t{tu}\t{tl}\t{wik}\t{total}\n")
EOF

############################################
# STEP 3: Assign contigs
############################################
echo "== Step 3: assigning contigs =="

python3 <<EOF
threshold = $THRESHOLD
min_kmers = $MIN_KMERS

assignments = {}

with open("$OUT/results/contig_counts.tsv") as f:
    next(f)
    for line in f:
        contig,ab,tu,tl,wik,total = line.strip().split()
        ab,tu,tl,wik,total = map(int,[ab,tu,tl,wik,total])

        if total < min_kmers:
            assignments[contig] = "ambiguous"
            continue

        d = {"AB":ab,"TU":tu,"TL":tl,"WIK":wik}
        best = max(d, key=d.get)
        frac = d[best] / total

        if frac >= threshold:
            assignments[contig] = best
        else:
            assignments[contig] = "ambiguous"

with open("$OUT/results/contig_assignments.tsv","w") as out:
    for c,a in assignments.items():
        out.write(f"{c}\t{a}\n")
EOF

############################################
# STEP 4: Split FASTA
############################################
echo "== Step 4: splitting FASTA =="

python3 <<EOF
assign = dict(line.strip().split() for line in open("$OUT/results/contig_assignments.tsv"))

files = {b:open(f"$OUT/fasta/{b}.fasta","w") for b in ["AB","TU","TL","WIK","ambiguous"]}

with open("$ASSEMBLY") as f:
    name=None; seq=[]
    for line in f:
        if line.startswith(">"):
            if name:
                b = assign.get(name,"ambiguous")
                files[b].write(f">{name}\n{''.join(seq)}\n")
            name=line[1:].split()[0]
            seq=[]
        else:
            seq.append(line.strip())

    if name:
        b = assign.get(name,"ambiguous")
        files[b].write(f">{name}\n{''.join(seq)}\n")

for f in files.values():
    f.close()
EOF

############################################
# STEP 5: Summary
############################################
echo "== Step 5: summary =="

cut -f2 $OUT/results/contig_assignments.tsv | sort | uniq -c \
| awk '{print $2"\t"$1}' > $OUT/results/summary.tsv


############################################
# STEP 6: Prepare Merqury database
############################################
echo "== Step 6: Preparing Merqury database =="

# Combine all grandparents k-mers (trusted set)
meryl union-sum output $OUT/merqury/all_parents.meryl \
    AB.meryl TU.meryl TL.meryl WIK.meryl

############################################
# STEP 7: Run Merqury per bin
############################################
echo "== Step 7: Running Merqury per strain bin =="

mkdir -p $OUT/merqury/results

for S in AB TU TL WIK ambiguous; do

    FASTA=$OUT/fasta/${S}.fasta

    # skip empty files
    if [ ! -s "$FASTA" ]; then
        echo "Skipping $S (empty)"
        continue
    fi

    echo "Running Merqury for $S"

    merqury \
        $OUT/merqury/all_parents.meryl \
        $FASTA \
        $OUT/merqury/results/${S}

done

echo "Pipeline complete!"
