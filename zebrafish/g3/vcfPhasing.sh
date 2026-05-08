#!/bin/bash

# Path to Beagle JAR file
beagle_jar="path/to/beagle.r1399.jar"

# Input VCF file
input_vcf="path/to/your/input.vcf"

# Output phased VCF file
output_vcf="path/to/your/output_phased.vcf"

# Run Beagle for phasing
java -Xmx4g -jar $beagle_jar \
  gt=$input_vcf \
  out=$output_vcf \
  impute=true

echo "Phasing completed. Phased VCF saved to $output_vcf"

