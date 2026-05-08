#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=winnowmap

#load the modules

module load winnowmap
module load samtools
module load meryl

#move to the right dir
cd  /data/okendojo/zebrafish/data/AB/asm/ont_ul_asm/asm_graph_map

graph="/data/okendojo/zebrafish/data/AB/asm/ont_ul_asm/5-untip/assembly.homopolymer-compressed.fasta"
ont="/data/okendojo/zebrafish/data/AB/batches_ont/ont_ul_hpc.fasta"


winnowmap -t 24  -cx map-ont ${graph} ${ont} | sed s/de:f://g  | awk -F "\t" '{ if ($12 >= 20 && $4-$3 > 5000 && 1-$21 >= 0.9) { if (match($5, "-")) print $1"\t"$2"\t"$3"\t"$4"\t+\t<"$6"\t"$7"\t"$7-$9"\t"$7-$8"\t"$10"\t"$11"\t"$12"\t"$13"\t"$15"\tdv:f:"$21"\tid:f:"1-$21; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t>"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$15"\tdv:f:"$21"\tid:f:"1-$21 }}' > winnowmap_ont_map.gaf 
