#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=750g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH	--gres=lscratch:1000
#SBATCH --job-name=mhc_cactus


#Load the modules
module load cactus

cd  /data/okendojo/zebrafish/data/g3/sex_project/mhcgene

rm -rf danio_mhc_pg

cactus-pangenome /lscratch/${SLURM_JOB_ID}/jobStore  /vf/users/okendojo/zebrafish/data/g3/sex_project/mhcgene/sequence.txt  --workDir /data/okendojo/zebrafish/data/g3/sex_project/mhcgene/workdir --outDir danio_mhc_pg --reference GRCz12tu --outName danio_mhc_pg  --snarlStats full --vcf filter --viz --giraffe --gfa  --gbz --odgi
