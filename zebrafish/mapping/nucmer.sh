#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=nucmer

#load the modules
module load mummer
module load syri
module load plotsr

cd /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/pur

nucmer --prefix=NHGRI_Fish6_GRCz11 NHGRI_Fish6_cons.fasta GRCz11.fasta
 

delta-filter -i 90  NHGRI_Fish6_GRCz11.delta > filtered.delta     # Remove small and lower quality alignments

show-coords -rcl filtered.delta > filtered.coords      # Convert alignment information to a .TSV format as required by SyRI

#syri -c GRCz11_fish11.filtered.coords --prefix syri_GRCz11_fish11 -d GRCz11_fish11.filtered.delta -r chr1_25.fasta -q fish11.fasta

#plotsr syri_GRCz11_fish11.out -o nucsynt.pdf  GRCz11_fish11 chr1_25.fasta fish11.fasta -H 8 -W 10
