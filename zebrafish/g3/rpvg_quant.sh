#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=200g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=rpvg_quant

cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/vcf_data/completeData/T12_quant

rpvg -t 24 -g "03.pantranscriptome/T12_quant.spliced.xg" -p  "03.pantranscriptome/T12_quant.haplotx.gbwt" -f "03.pantranscriptome/T12_quant.txorigin.tsv" -a "04.multimapping/T12_quant.aligned.gamp" -o time12_hst_quantification --inference-model haplotype-transcripts

