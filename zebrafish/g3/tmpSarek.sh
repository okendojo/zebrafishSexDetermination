#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=210g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=sarek_variantcall

#load modules
module add nextflow
module add singularity


cd  /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis


nextflow run nf-core/sarek --input 'samplesheet.csv'  --genome null   --skip_tools haplotypecaller_filter --skip_tools haplotyper_filter --step 'markduplicates'  --igenomes_ignore  --fasta '/data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.dna.primary_assembly.fa'   --fasta_fai '/data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.dna.primary_assembly.fa.fai' --outdir 'rna_timepoint_variants_GATK'   --tools 'haplotypecaller'   --trim_fastq true   --concatenate_vcfs false   --aligner 'dragmap' --known_snps '/data/okendojo/zebrafish/data/g3/illumina/AB_F_CB/qc/dv_WGS/dbSNP.vcf.gz' --skip_tools baserecalibrator --known_snps_tbi '/data/okendojo/zebrafish/data/g3/illumina/AB_F_CB/qc/dv_WGS/dbSNP.vcf.gz.tbi'  -resume  -profile biowulf
