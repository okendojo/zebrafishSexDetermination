#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=240g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=LASTZ_ALIGNMENT

#load the modules

module load winnowmap
module add delly
module load samtools
module add lastz

cd /data/okendojo/mate_data

#lastz chr9.fa  ich_contigs.fasta --notransition --step=20 --nogapped --identity=100 --seed=14of22 --output=ich_alignment_output.maf --format=maf 


lastz chr9.fa  ich_contigs.fasta --output=ich_Test_alignment_output.maf --format=maf

#ich_contigs.fasta[multiple]
