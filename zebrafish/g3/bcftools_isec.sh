#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=120g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=bcftools_isec

#load module
module load bcftools

cd /data/okendojo/zebrafish/data/g3/illumina/deepvariantVCFs

#!/bin/bash

# Paths to the VCF files (replace these with your actual file paths)
AB="AB.vcf.gz"
WIK="WIK.vcf.gz"
TL="TL.vcf.gz"
TU="TU.vcf.gz"

# Output directory
output_dir="variant_comparison_output"
mkdir -p "$output_dir"

# Compare VCF files and extract unique variants
bcftools isec -n=4 -c all -Oz \
  -o "$output_dir/unique_variants_all.vcf.gz" \
  "$AB" "$WIK" "$TL" "$TU"

# Extract unique variants from each file
for i in {1..4}; do
  bcftools view -Oz -o "$output_dir/unique_variants_file${i}.vcf.gz" \
    "$output_dir/0000.vcf.gz" "$output_dir/000${i}.vcf.gz"
done

echo "Comparison completed. Unique variants for each file are in $output_dir."


#The bcftools isec command is used to compare the VCF files and extract variants that are unique to each file. The -n=4 option specifies that all four VCF files are being compared, and -c all ensures that variants present in all files are excluded from the output. The bcftools view command is then used to extract the unique variants from each file individually.
