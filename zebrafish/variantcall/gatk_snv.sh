#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=64g
#SBATCH --ntasks-per-core=1
#SBATCH --time=96:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov

# Workflow is adopted from: https://gatk.broadinstitute.org/hc/en-us/articles/4407897446939--How-to-Run-germline-single-sample-short-variant-discovery-in-DRAGEN-mode#article-comments
#Load the required module
module load GATK
module load dragmap
#module load singularity
#move to the right dir
cd  /data/okendojo/paradisfishProject/snv

#DRAGMAP - alignment step; Build new hash table of a reference file

REF=/data/okendojo/paradisfishProject/snv/reference/macOpe2Assembly.fasta
OUTDIR=/data/okendojo/paradisfishProject/snv/reference
FQ1=/data/okendojo/paradisfishProject/snv/ERR3332352_1.fastq.gz
FQ2=/data/okendojo/paradisfishProject/snv/ERR3332352_2.fastq.gz
STRTABLE=/data/okendojo/paradisfishProject/snv/str_table.tsv


#dragen-os --build-hash-table true --ht-reference ${REF} --output-directory ${OUTDIR}

#Map reads to the reference with DRAGMAP
dragen-os -r /data/okendojo/paradisfishProject/snv/dragenRef -1 ERR3332352_1.fastq.gz -2 ERR3332352_2.fastq.gz --num-threads=64 --ht-num-threads=64 --ht-mem-limit=128GB  --output-directory /data/okendojo/paradisfishProject/snv/mappedData --output-file-prefix macOpe2 

#Data pre-processing; Create STR Table File for the Reference
#gatk ComposeSTRTableFile -R reference.fasta -O str_table.tsv

#Build the DRAGEN STR Model from the Aligned Reads
#gatk CalibrateDragstrModel -R ${REF} -I macOpe2.sam -str ${STRTABLE} -O dragstr_model.txt

#Variant Calling and Filtering
#gatk HaplotypeCaller -R ${REF} -I macOpe2.sam -O macOpe2_snvs.vcf --dragen-mode true --dragstr-params-path dragstr_model.txt

#Hard Filter Variants
#gatk VariantFiltration -V macOpe2_snvs.vcf --filter-expression "QUAL < 10.4139" --filter-name "DRAGENHardQUAL" -O macOpe2_snv_filtered.vcf
















