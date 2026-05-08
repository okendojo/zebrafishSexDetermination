#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=800g
#SBATCH --gres=lscratch:500
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=F11_hifiasm

#source /data/$USER/conda/etc/profile.d/conda.sh && source /data/$USER/conda/etc/profile.d/mamba.sh

#mamba activate hifiasm #verkko and hifiasm are installed together

#module load
module load hifiasm/0.19.5

export TMPDIR=/lscratch/$SLURM_JOB_ID

cd /data/okendojo/zebrafish/data/fish11/hifiasm


#Run the assembler 
hifiasm -o asm.out -t 64 -l 0 --ul /data/okendojo/zebrafish/data/fish11/ont/fastq_pass/concatenatedONT/zfish11ONT.fastq.gz  /data/okendojo/zebrafish/data/fish11/hifi/*.fastq.gz	
