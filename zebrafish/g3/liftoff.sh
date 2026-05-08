#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=liftoff


#Load the modules
#source /data/$USER/conda/etc/profile.d/conda.sh && source /data/$USER/conda/etc/profile.d/mamba.sh
source myconda
#conda activate liftoff
#module load python/3.9
module load minimap2

cd /data/okendojo/paradisfishProject/genoComparison/splendent_liftoff


#liftoff  /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/asm_polished/NHGRI_Fish6_cons.fasta  /data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.dna.primary_assembly.fa -g /data/okendojo/zebrafish/refGenome/Danio_rerio.GRCz11.111.gtf  -o GRCz12tu.gtf  -u unmapped_features2.txt  -copies -polish -m minimap2 -dir liftof_tmp2 -chroms /data/okendojo/zebrafish/data/fish6/asm/t2t/liftoff/ensemble/chrs.txt

liftoff macOpe1_asm.fasta GCF_900634795.4_fBetSpl5.4_genomic.fna -g GCF_900634795.4_fBetSpl5.4_genomic.gtf -o macOpe2.gtf -u unmapped_features2.txt  -copies -p 24 -polish -m minimap2 -flank 0.8 -dir liftof_tmp2 -chroms chrs.txt 
