#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=230g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH	--gres=lscratch:1000
#SBATCH --job-name=cb_ekk_pg


#Load the modules
module load cactus
ml singularity
ml docker

cd /data/okendojo/zebrafish/data/g3/sex_project/zebrafishpg 

cactus-pangenome /lscratch/${SLURM_JOB_ID}/jobStore /data/okendojo/zebrafish/data/g3/sex_project/asm/fa/sequencefiles.txt  --workDir tmpwkdir2 --outDir cb_ekk_pg --reference GRCz12tu --outName drerio_pg --snarlStats full --vcf --viz --giraffe  --gfa --chrom-vg --chrom-og --odgi --cleanWorkDir always --haplo --gbz --vcfwave
