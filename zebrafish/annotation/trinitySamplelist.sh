#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=500g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=trinity_asm


cd /data/okendojo/zebrafish/data/g3/rna_sequences/trimmed

# Load the module
module load trinity/2.14.0

mkdir trinityasm_zfish

Trinity --seqType fq --max_memory 450G --full_cleanup  --samples_file sample.txt  --monitoring --monitor_sec 30 --CPU 32 --output trinityasm_zfish 


