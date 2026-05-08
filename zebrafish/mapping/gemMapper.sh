#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=T2T_GEM_MAPPING_sam

#load the modules

module load winnowmap
module add delly
module load samtools

#paths to the genomic data
R1="/data/Zebrafish_T2T/fish11/illum/R1.fastq.gz"
R2="/data/Zebrafish_T2T/fish11/illum/R2.fastq.gz"
T2T="/data/Zebrafish_T2T/fish11/resolved_assembly/fish11_merfin_polished.fasta"


#move to the right dir
WKDIR="/data/okendojo/zebrafish/data/centromere_analysis/satelliteRepFinder/fish11/illumina_map"

cd ${WKDIR}

#Index the reference genome
gem-indexer -i ${T2T} -o fish11_T2T

#Paired-end mapping
gem-mapper -I fish11_T2T.gem -1 $R1 -2 $R2 -z --report-file=analysis_report_sam -p  -t 24 -o fish11_illumina_pe.sam 


samtools view -S -b fish11_illumina_pe.sam -o fish11_illumina_pe.bam
samtools sort fish11_illumina_pe.bam -o fish11_illumina_pe_sorted.bam
samtools index fish11_illumina_pe_sorted.bam


rm fish11_illumina_pe.sam



