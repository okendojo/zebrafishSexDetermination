#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=F11_resolvedGraph

#load module
module add minimap2
module add bedtools
module add samtools
module add graphaligner
module add seqtk
module add python
module add snakemake
module add R
#move to the dir
cd  /data/okendojo/zebrafish/data/fish11/sg_sandbox

./src/resolution/use_ont.sh tangle-resolution/asm_graphs/ont_hpc.fasta tangle-resolution/asm_graphs/assembly.homopolymer-compressed.gfa contigs.txt resolvedGraphs
