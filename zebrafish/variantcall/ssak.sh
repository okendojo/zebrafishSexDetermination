#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=Temp_variantCall

#Load the required modules
module purge
module load nextflow
module load singularity


#move to the directory containing the data
cd /data/okendojo/zebrafish/data/g3


CONFIG=/data/okendojo/slugProject/annot/scrna/nextflow.config
FASTA='/data/okendojo/zebrafish/refGenome/GRCz11_genomic.fasta'
FAI='/data/okendojo/zebrafish/refGenome/GRCz11_genomic.fasta.fai'
KNOWNSNPS='/data/okendojo/zebrafish/data/g3/known_snps/GCA_000002035.4_current_ids.vcf.gz'
SNPTBI='/data/okendojo/zebrafish/data/g3/known_snps/GCA_000002035.4_current_ids.vcf.gz.tbi' 


mkdir g3_temp2
#--skip_tools baserecalibrator 
nextflow run sarek --input samplesheet.csv -resume --outdir g3_temp2  --tools 'freebayes,manta,cnvkit,vep' --dbsnp '/data/okendojo/zebrafish/data/g3/known_snps/GCA_000002035.4_current_ids.vcf.gz' --dbsnp_tbi '/data/okendojo/zebrafish/data/g3/known_snps/GCA_000002035.4_current_ids.vcf.gz.tbi' --genome GRCz10 -profile singularity -c ${CONFIG} --known_snps $KNOWNSNPS --known_snps_tbi $SNPTBI 
