#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=230g
#SBATCH --ntasks-per-core=1
#SBATCH --time=140:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=cb_verkko


#Load the modules

module load verkko/2.2.1
module load python/3.10
module load seqkit

#Run the reads extraction first
cd  /data/okendojo/zebrafish/data/g3/sex_project/asm

#Run verkko

verkko -d asm_cb --hifi  /data/okendojo/zebrafish/data/g3/sex_project/fastq_files/coochbehar/*.gz --no-nano --grid  --lay-run 16 100 96  --ali-run 16 100 96 
