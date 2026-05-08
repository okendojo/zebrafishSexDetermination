#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=232g
#SBATCH --ntasks-per-core=1
#SBATCH --time=170:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=ZFish_repmodeler

#Load the required modules
module add repeatmasker
module add repeatmodeler

cd /data/okendojo/zebrafish/data/fish11/annotation

# build new RepeatModeler BLAST database
BuildDatabase -name zebrafish -engine ncbi fish11_merfin_polished.fasta

# now run RepeatModeler with 32 cores (we have alot of resourcee), you may have to scale it according to your resources
RepeatModeler -pa 32 -engine ncbi -LTRStruct -recoverDir ./ -database zebrafish 2>&1 | tee 00_repeatmodeler.log
