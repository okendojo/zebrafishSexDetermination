#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=trio_verkko


#Load the modules
module load verkko/2.2


#Run the reads extraction first

cd  /data/okendojo/zebrafish/data/g3/assembly

#Run verkko

#verkko -d asm_trio4gps --hifi  /data/okendojo/zebrafish/data/g3/pacBio/hifi.fastq.gz --nano /data/okendojo/zebrafish/data/g3/ontData/combined/deduped_ont.fasta --hap-kmers /data/okendojo/zebrafish/data/ab_asm/binning/dfraw/hapmers/grandparents/hapmers/AB_TU.only.meryl /data/okendojo/zebrafish/data/ab_asm/binning/dfraw/hapmers/grandparents/hapmers/WIK_TL.only.meryl trio --slurm 


verkko -d AB_TU --hifi  /data/okendojo/zebrafish/data/g3/binning/tu_ab/pacbio_classified/*.fastq.gz --nano /data/okendojo/zebrafish/data/g3/binning/tu_ab/ont_classified/ont.fa --hap-kmers /data/okendojo/zebrafish/data/ab_asm/binning/dfraw/hapmers/finalcompressedhapmer/AB.k21.only.meryl  /data/okendojo/zebrafish/data/ab_asm/binning/dfraw/hapmers/finalcompressedhapmer/TU.k21.only.meryl trio --slurm
