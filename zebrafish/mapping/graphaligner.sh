#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=ab_grapaligner


#Load the modules
module load mbg
module load minimap2
module load graphaligner
module load mashmap
module load winnowmap
module load singularity
ml seqtk seqkit
ml seqtk

cd /data/okendojo/zebrafish/data/g3/output/meryl/bins/abtu/ont

seqtk hpc ont.fasta.gz > ont_hpc.fasta.gz

cd  /data/okendojo/zebrafish/data/g3/assembly/asm_abtu

set -e

if [ 5000 -gt 0 ]; then
   memwindow="--seeds-mxm-windowsize 5000"
fi

ont="/data/okendojo/zebrafish/data/g3/output/meryl/bins/abtu/ont/ont_hpc.fasta.gz"
graph="assembly.homopolymer-compressed.gfa"

GraphAligner -t 24 -g ${graph} -f ${ont} -a ab_graphaliner.gaf \
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
