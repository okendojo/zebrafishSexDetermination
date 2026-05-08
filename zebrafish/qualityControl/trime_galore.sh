#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=RNA_fastqTrim

# Load Trim Galore! module
module load trimgalore
module load multiqc
module load fastqc

cd /data/okendojo/zebrafish/data/g3/rna_sequences

# Define the path to the FastQ files directory
input_dir="/data/okendojo/zebrafish/data/g3/rna_sequences/untrimmed"

# Define the output directory
output_dir="/data/okendojo/zebrafish/data/g3/rna_sequences/trimmed"

# Loop through all FastQ files in the input directory
for fastq_file in "$input_dir"/*.fastq.gz; do
  # Extract the filename without the extension
  filename=$(basename "$fastq_file" .fastq.gz)

  # Run Trim Galore! on the FastQ file
  trim_galore --quality 20  --fastqc_args "--threads 24 --min_length 25" --cores 24 --output_dir "$output_dir" "$fastq_file"

done

cd $output_dir

multiqc --title 'F2_RNAseq quality control' . 
