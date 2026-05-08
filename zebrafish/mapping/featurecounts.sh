#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=featurecounts

#load the modules
module load subread

cd /data/okendojo/zebrafish/data/g3/assembly/asm_mapping

GTF="/data/okendojo/zebrafish/refGenome/three_gen_asm.gtf"

featureCounts -T 24 --primary -t exon -g gene_id -a $GTF -o g3_featurecounts.txt GRCZ11_G3_mapping_sorted.bam 

