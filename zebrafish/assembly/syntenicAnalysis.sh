
#Align Genomes Using Minimap2
minimap2 -x asm5 -t 4 genome1.fa genome2.fa > genome1_vs_genome2.paf

#Convert PAF to MAF (Multiple Alignment Format) for Synteny Analysis
paftools.js view -f maf genome1_vs_genome2.paf > genome1_vs_genome2.maf

#Visualize Synteny Blocks
conda install -c bioconda synteny_plot
synteny_plot --maf genome1_vs_genome2.maf --output synteny_plot.png

