#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=dellyVCFCall

#load the modules

module load winnowmap
module add delly
module load samtools

#paths to the genomic data
hifi="/data/okendojo/zebrafish/data/fish11/hifi/*.fastq.gz"
ref="/data/Zebrafish_T2T/fish11/resolved_assembly/fish11_merfin_polished.fasta"
grzz11="/data/okendojo/zebrafish/data/fish11/gggenome/original/ref.fasta"
#move to the right dir
WKDIR="/data/okendojo/zebrafish/data/centromere_analysis/satelliteRepFinder/fish11"

cd ${WKDIR}


#mapping PacBio hifi data

#winnowmap -t 24 -ax map-pb ${ref} ${hifi} -o fish11_pbMapping.sam


echo "Now running samtools.........."

#samtools view -S -b fish11_pbMapping.sam -o fish11_pbMapping.bam
#samtools sort fish11_pbMapping.bam -o fish11_pbMapping_sorted.bam
#samtools index fish11_pbMapping_sorted.bam

echo "Now running delly!"

delly lr -y pb -o delly_fish11_grcz11.bcf -q 20 -t ALL -g ${ref} fish11_pbMapping_sorted.bam 


