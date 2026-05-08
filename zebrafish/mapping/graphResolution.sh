#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=graphresolution

#load module
module add minimap2
module add bedtools
module add samtools
module add graphaligner
module add seqtk
module add python
module add snakemake
module add R
#move to the dir
cd  /data/okendojo/zebrafish/data/fish11/sg_sandbox

./src/canu_launch/master_trim.sh /data/okendojo/zebrafish/data/fish11/sg_sandbox/src/pipe/config.yaml  tmpGraph 150000 /data/okendojo/zebrafish/data/fish11/hifi/hifi.fasta
