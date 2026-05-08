#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=90:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=starMap

#load modules
module add nextflow
module add singularity
module add STAR

cd /data/okendojo/zebrafish/data/g3/rna_sequences/untrimmed

# define variables

genome_index="/data/okendojo/zebrafish/data/g3/rna_sequences/gene_quantification/rnaseq_result/star_index"


# Output directory
output_dir="/data/okendojo/zebrafish/data/g3/eQTL/BAMs"

# Array of sample names
samples=("22CMMHLT3_19570451_S60_L003" "22CMMHLT3_19570453_S61_L003" "22CMMHLT3_19570455_S62_L003" "22CMMHLT3_19570457_S63_L003" "22CMMHLT3_19570459_S64_L003" "22CMMHLT3_19570461_S65_L003" "22CMMHLT3_19570463_S66_L003" "22CMMHLT3_19570465_S67_L003" "22CMMHLT3_19570467_S68_L003" "22CMMHLT3_19570469_S69_L003" "22CMMHLT3_19570471_S70_L003")

# Loop through samples and map paired-end RNA-seq reads
for sample in "${samples[@]}"; do
    # Paths to paired-end reads for the current sample
    reads1="/data/okendojo/zebrafish/data/g3/rna_sequences/untrimmed/${sample}_R1_001.fastq.gz"
    reads2="/data/okendojo/zebrafish/data/g3/rna_sequences/untrimmed/${sample}_R2_001.fastq.gz"

    # Output directory for the current sample
    sample_output_dir="$output_dir/$sample"

    # Create the output directory
    mkdir -p "$sample_output_dir"

    # Run STAR aligner
    STAR \
        --runThreadN 32 \
	--readFilesCommand zcat \
	--twopassMode Basic \
        --genomeDir $genome_index \
        --readFilesIn $reads1 $reads2 \
        --outFileNamePrefix "$sample_output_dir/" \
        --outSAMtype BAM SortedByCoordinate

    echo "Mapping for $sample complete."
done

echo "Mapping of paired-end RNA-seq reads complete."

