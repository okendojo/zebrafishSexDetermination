#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=1000g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=only_meryl_lookup

# Load meryl module
module add meryl

#move the directory containing the data
cd /data/okendojo/zebrafish/data/g3/eQTL/dnaVariants/deepvariantStrainSpecificvariants/variants/secondAnalysis/vcf_data/completeData


meryl-lookup -existence -sequence /data/okendojo/zebrafish/data/g3/rna_sequences/untrimmed/t24_hpc_R1.fasta  -mers /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/AB.only.meryl/ /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/TU.only.meryl/ /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/TL.only.meryl/ /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/WIK.only.meryl/ -output t24_forward_contigKemrs.txt 


meryl-lookup -existence -sequence /data/okendojo/zebrafish/data/g3/rna_sequences/untrimmed/t24_hpc_R2.fasta -mers /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/AB.only.meryl/ /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/TU.only.meryl/ /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/TL.only.meryl/ /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/WIK.only.meryl/ -output t24_reverse_contigKemrs.txt

#meryl-lookup -existence -sequence rna_seq.fq.gz -mers $dataPath/TU.only.meryl $dataPath/WIK.only.meryl $dataPath/AB.only.meryl $dataPath/TL.only.meryl > RNA-Seq_only.hapmer_count.txt
