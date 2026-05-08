#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=rDNA_array

#Load modules
module load samtools
module load blast
module load bedtools

cd /data/okendojo/zebrafish/data/zfish_rDNA

# Input files
genome="/data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/asm_polished/NHGRI_Fish6_cons.fasta"
rDNA_query="/data/okendojo/zebrafish/data/zfish_rDNA/rDNA_sequence.fasta"  # Provide rDNA sequences (e.g. 18S, 5.8S, 28S, 5S)

# Output directories
mkdir -p results/blast results/bed

# Step 1: Index the genome
samtools faidx $genome

# Step 2: Create BLAST database
makeblastdb -in $genome -dbtype nucl -out results/blast/genome_db

# Step 3: Run BLAST to find rDNA hits
blastn -query $rDNA_query -db results/blast/genome_db \
  -out results/blast/rDNA_hits.tsv -outfmt 6 -evalue 1e-10

# Step 4: Convert BLAST hits to BED
awk '{start=($9<$10)?$9:$10; end=($9>$10)?$9:$10; print $2"\t"start"\t"end}' \
  results/blast/rDNA_hits.tsv | sort -k1,1 -k2,2n > results/bed/rDNA_hits.bed

# Step 5: Merge nearby hits to cluster rDNA arrays
bedtools merge -i results/bed/rDNA_hits.bed -d 1000 > results/bed/rDNA_clusters.bed

# Step 6: Optional - Visualize in IGV or convert to Circos format using custom scripts
# IGV: Load genome.fasta and results/bed/rDNA_clusters.bed
# Circos: Use results/bed/rDNA_clusters.bed to create karyotype or link tracks

