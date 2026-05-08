#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=230g
#SBATCH --ntasks-per-core=1
#SBATCH --time=140:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=G3_verkko


#Load the modules

module load verkko/2.2.1
module load python/3.10
module load seqkit

#Run the reads extraction first
cd  /data/okendojo/zebrafish/data/g3/assembly

#Run verkko
#trio binning
#verkko -d asm_120325 --hifi  /data/okendojo/zebrafish/data/g3/pacBio/*.gz --nano  /data/okendojo/zebrafish/data/g3/ontData/ont_uniq.fastq.gz --hap-kmers /data/okendojo/zebrafish/data/g3/merylDBs/f1gen/hapmers/maternal.k30.hapmer.meryl /data/okendojo/zebrafish/data/g3/merylDBs/f1gen/hapmers/paternal.k30.hapmer.meryl trio --grid  --lay-run 16 100 96  --ali-run 16 100 96 


verkko -d asm_jan2026 --hifi  /data/okendojo/zebrafish/data/g3/pacBio/*.gz  --nano  /data/okendojo/zebrafish/data/g3/ontData/ONT.fastq.gz --grid  --lay-run 16 100 96  --ali-run 16 100 96 
