#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=zip

cd /data/okendojo/zebrafish/data/g3

cp -r iso-seq/ /data/okendojo/datashare/quanrong_model/long_reads/

cd /data/okendojo/zebrafish/data/g3/rna_sequences/untrimmed/TX

cp * /data/okendojo/datashare/quanrong_model/illumina/

cd /data/okendojo/datashare/ 

tar -cvzf quanrong_model.tar.gz quanrong_model

rm quanrong_model 
