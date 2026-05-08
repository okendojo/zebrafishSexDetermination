#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=f6_graphAlign


#Load modulels
module load seqkit
module load seqtk
module load graphaligner

#move to the working dir

cd /data/okendojo/zebrafish/data/fish6/asm/alignmentGraph


#rm ont_reads.fasta

#Run graphaligner
set -e

if [ 5000 -gt 0 ]; then
   memwindow="--seeds-mxm-windowsize 5000"
fi

graph="/data/okendojo/zebrafish/data/fish6/asm/verkko_asm_fish6/assembly.homopolymer-compressed.gfa"
ont="/data/okendojo/zebrafish/data/fish6/ontData/concatenated/ont_hpc.fasta"

GraphAligner -t 24 -g ${graph} -f ${ont} -a ont_Graphaligned.gaf \
  $memwindow --seeds-mxm-length 30 \
  --seeds-mem-count 10000 \
  --bandwidth 15 \
  --multimap-score-fraction 0.99 \
  --precise-clipping 0.85 \
  --min-alignment-score 5000 \
  --hpc-collapse-reads \
  --discard-cigar \
  --clip-ambiguous-ends 100 \
  --overlap-incompatible-cutoff 0.15 \
  --max-trace-count 5 \
  --mem-index-no-wavelet-tree

