#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --gres=lscratch:500
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=zfishRepeatIDs


#Load the required modules
module load maker
module load augustus
module load busco
module load snap
module load bedtools
module load repeatmodeler
#module load repeatmasker
module load genometools
module load seqkit

cd /data/okendojo/zebrafish/data/fish11/annotation

#1. De Novo Repeat Identification
# build new RepeatModeler BLAST database
BuildDatabase -name zebrafish -engine ncbi /data/okendojo/zebrafish/data/shawmasmres/fish11_hifi_ul_verkko.fasta

# now run RepeatModeler with 32 cores (we have alot of resourcee), you may have to scale it according to your resources
RepeatModeler -pa 32 -engine ncbi -LTRStruct -database zebrafish 2>&1 | tee 00_repeatmodeler.log

