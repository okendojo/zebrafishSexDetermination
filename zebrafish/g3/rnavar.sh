#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=200g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=rnavar

#load modules
module add nextflow
module add singularity


cd /data/okendojo/zebrafish/data/g3/F2_variants

fasta="/data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/asm_polished/NHGRI_Fish6_cons.fasta"

nextflow run nf-core/rnavar --input 'samplesheet.csv'  --genome null  --igenomes_ignore  -c "nextflow.config"  --fasta ${fasta} --skip_baserecalibration false --outdir 'f2_rna_variants' --dbsnp '/data/okendojo/zebrafish/data/g3/F2_variants/dbsnps/db1.vcf.gz' --dbsnp_tbi '/data/okendojo/zebrafish/data/g3/F2_variants/dbsnps/db1.vcf.gz.tbi' --gtf '/data/okendojo/zebrafish/data/g3/F2_variants/fish6_fixedAgatfixed.gtf' --star_index 'grcz12_star_index'  --skip_baserecalibration false  --skip_intervallisttools false  --skip_variantfiltration false  --skip_variantannotation true  --skip_multiqc false -resume --remove_duplicates true  -profile biowulf

#rm -rf work .nextflow*
