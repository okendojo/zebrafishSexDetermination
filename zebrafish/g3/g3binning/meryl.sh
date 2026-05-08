#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=360g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=maryl


#Load the modules
ml meryl

cd /data/okendojo/zebrafish/data/g3/merylDBs/uncompressed

# Count k-mers for each founder
#meryl count k=31 /data/okendojo/zebrafish/data/g3/illumina/AB_F_CB/HFVYJDSX5_19412535_S54_L002_R1_001.fastq.gz /data/okendojo/zebrafish/data/g3/illumina/AB_F_CB/HFVYJDSX5_19412535_S54_L002_R2_001.fastq.gz  output AB.meryl
#meryl count k=31 /data/okendojo/zebrafish/data/g3/illumina/TL_F_CB/HFVYJDSX5_19412527_S51_L002_R1_001.fastq.gz /data/okendojo/zebrafish/data/g3/illumina/TL_F_CB/HFVYJDSX5_19412527_S51_L002_R2_001.fastq.gz output TL.meryl
#meryl count k=31 /data/okendojo/zebrafish/data/g3/illumina/TU_M_CB/HFVYJDSX5_19412541_S50_L002_R1_001.fastq.gz /data/okendojo/zebrafish/data/g3/illumina/TU_M_CB/HFVYJDSX5_19412541_S50_L002_R2_001.fastq.gz output TU.meryl
#meryl count k=31 /data/okendojo/zebrafish/data/g3/illumina/WIK_M_CB/HFVYJDSX5_19412539_S55_L002_R1_001.fastq.gz /data/okendojo/zebrafish/data/g3/illumina/WIK_M_CB/HFVYJDSX5_19412539_S55_L002_R2_001.fastq.gz output WIK.meryl

# Make "unique" k-mer sets per founder (present in founder, absent in the other three)
#meryl difference AB.meryl TL.meryl TU.meryl WIK.meryl output AB.unique.meryl
#meryl difference TL.meryl AB.meryl TU.meryl WIK.meryl output TL.unique.meryl
#meryl difference TU.meryl AB.meryl TL.meryl WIK.meryl output TU.unique.meryl
#meryl difference WIK.meryl AB.meryl TL.meryl TU.meryl output WIK.unique.meryl

#scan the assembly

#meryl-lookup -existence -mers AB.unique.meryl -sequence /data/okendojo/zebrafish/data/g3/assembly/asm_270625/assembly.fasta > AB.hits.txt
#meryl-lookup -existence -mers TL.unique.meryl -sequence /data/okendojo/zebrafish/data/g3/assembly/asm_270625/assembly.fasta > TL.hits.txt
#meryl-lookup -existence -mers TU.unique.meryl -sequence /data/okendojo/zebrafish/data/g3/assembly/asm_270625/assembly.fasta > TU.hits.txt
#meryl-lookup -existence -mers WIK.unique.meryl -sequence /data/okendojo/zebrafish/data/g3/assembly/asm_270625/assembly.fasta > WIK.hits.txt

meryl-lookup -existence -sequence /data/okendojo/zebrafish/data/g3/assembly/asm_270625/assembly.fasta -mers AB.unique.meryl TU.unique.meryl TL.unique.meryl WIK.unique.meryl > f2_genKmers.txt
