#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=meryl

#load the modules

module load winnowmap
module load samtools
module load meryl

#Run the pipeline

cd /data/Zebrafish_T2T/fish11/illum

meryl count k=21  *_001.fastq.gz output illumina.meryl

mv illumina.meryl /data/okendojo/zebrafish/data/fish11/polish/mapping_polished/hybridMeryl_illumHiFi


cd /data/okendojo/zebrafish/data/fish11/polish/mapping_polished/hybridMeryl_illumHiFi

meryl count k=21 /data/okendojo/zebrafish/data/fish11/hifi/*.fastq.gz output hifi.meryl

meryl greater-than 1 illumina.meryl output illumina_gt1.meryl
meryl greater-than 1 hifi.meryl output hifi_gt1.meryl 

meryl union-sum illumina_gt1.meryl hifi_gt1.meryl output hifi_illum_hybrid.meryl

FASTA=/vf/users/okendojo/zebrafish/data/fish11/polish/mapping_polished/fasta/fish11_merfin_polished.fasta

$tools/merqury/_submit_merqury.sh hifi_illum_hybrid.meryl ${FASTA} hybrid_hifi_illum_kmer
