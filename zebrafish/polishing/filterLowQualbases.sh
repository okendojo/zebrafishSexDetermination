#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=baseFilter

#load the modules

module load winnowmap
module load samtools
module load meryl
module add minimap2
module add merfin
module add bcftools
module add racon
module add pbipa 
module add bedtools

cd /data/okendojo/zebrafish/data/fish11/polishing/racon
# Define input genome assembly in FASTA format
genome_assembly="fish11.fasta"

# Define a threshold for sequence quality (e.g., minimum base quality or length)
quality_threshold=20

# Define the output directory
output_dir="low_quality_regions"

# Create the output directory
mkdir -p $output_dir

# Step 1: Calculate sequence quality (e.g., base quality or length) across the genome
# In this example, we calculate base quality and consider bases below the threshold as low-quality.
bedtools nuc -fi $genome_assembly | \
  awk -v threshold=$quality_threshold '$4 < threshold {print $1 "\t" $2-1 "\t" $2}' \
  > $output_dir/low_quality_regions.bed

# Step 2: Optionally, you can annotate these low-quality regions using BEDTools intersect or other tools.
# For example, if you have information about genes or features, you can intersect the low-quality regions with gene annotations.

# Step 3: Generate a summary report or visualization of low-quality regions (optional).

# Done!
echo "Low-quality genome assembly regions cataloged in $output_dir."

