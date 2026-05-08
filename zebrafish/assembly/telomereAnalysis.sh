
#Run telomere analysis script (VGP) script. This will create bed file to be used in the subsequent sections

$telomere/pipeline/telomere/telomere_analysis.sh zebrafish_210xT2T 0.4 50 assembly.fasta

#run mash to determine rDNA coordinates
mash sketch rdna.fasta
mash screen rdna.fasta.msh assembly.fasta > rdna.screennodes.out


#serge prefered telomere analysis script
/data/korens/devel/utils/telomere/find assembly.fasta > telomere

java -cp /data/korens/devel/utils:. FindTelomereWindows telomere 99.9 | awk '{if ($NF > 0.5) print $2"\t"$4"\t"$5"\t"$3"\t"$NF}' | sed s/\>//g|bedtools merge -d -500 -i - -c 4 -o distinct | bedtools sort -i - | bedtools merge -i - > tmp
mv tmp telomere

#Then run the following
python /data/korens/devel/verkko-fixgaps/src/scripts/remove_nodes_add_telomere.py assembly.homopolymer-compressed.noseq.gfa assembly.fasta.fai rdna.screennodes.out telomere assembly.scfmap > assembly.homopolymer-compressed.nordna.noseq.gfa 2> tmp.err

#This is the latest run
python /data/korens/devel/verkko-fixgaps/src/scripts/remove_nodes_add_telomere.py -r rdna.screennodes.out -t telomeres.bed -g assembly.homopolymer-compressed.noseq.gfa -s assembly.scfmap -o assembly.homopolymer-compressed.nordna.noseq.gfa

winnowmap -t 24 -cx map-ont graph.hpc.fasta reads.hpc.fasta | sed s/de:f://g | awk -F "\t" '{ if ($12 >= 20 && $4-$3 > 5000 && 1-$21 >= 0.9) { if (match($5, "-")) print $1"\t"$2"\t"$3"\t"$4"\t+\t<"$6"\t"$7"\t"$7-$9"\t"$7-$8"\t"$10"\t"$11"\t"$12"\t"$13"\t"$15"\tdv:f:"$21"\tid:f:"1-$21; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t>"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$15"\tdv:f:"$21"\tid:f:"1-$21 }}'
