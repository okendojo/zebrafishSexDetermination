#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=32
#SBATCH --mem=1000g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=meryl

# Load meryl module
module add seqtk
module load meryl

cd /data/okendojo/zebrafish/data/g3/rna_sequences/untrimmed

#combine timepoints

cat 22CMMHLT3_19570467_S68_L003_R1_001.fastq.gz 22CMMHLT3_19570469_S69_L003_R1_001.fastq.gz 22CMMHLT3_19570471_S70_L003_R1_001.fastq.gz > t24_forward_reads_R1.fastq.gz

cat 22CMMHLT3_19570467_S68_L003_R2_001.fastq.gz 22CMMHLT3_19570469_S69_L003_R2_001.fastq.gz 22CMMHLT3_19570471_S70_L003_R2_001.fastq.gz  > t24_reverse_reads_R2.fastq.gz



seqtk seq -A t24_forward_reads_R1.fastq.gz > t24_forward_reads_R1.fasta

seqtk seq -A t24_reverse_reads_R2.fastq.gz > t24_reverse_reads_R2.fasta

seqtk hpc t24_forward_reads_R1.fasta > t24_forward_reads_hpc_R1.fasta
seqtk hpc t24_reverse_reads_R2.fasta > t24_reverse_reads_hpc_R2.fasta

rm t24_forward_reads_R1.fasta t24_reverse_reads_R2.fasta 


cat 22CMMHLT3_19570461_S65_L003_R1_001.fastq.gz 22CMMHLT3_19570463_S66_L003_R1_001.fastq.gz 22CMMHLT3_19570465_S67_L003_R1_001.fastq.gz > t12_forward_reads_R1.fastq.gz

cat 22CMMHLT3_19570461_S65_L003_R2_001.fastq.gz 22CMMHLT3_19570463_S66_L003_R2_001.fastq.gz 22CMMHLT3_19570465_S67_L003_R2_001.fastq.gz > t12_reverse_reads_R2.fastq.gz


seqtk seq -A t12_forward_reads_R1.fastq.gz > t12_forward_reads_R1.fasta

seqtk seq -A t12_reverse_reads_R2.fastq.gz > t12_reverse_reads_R2.fasta

seqtk hpc t12_forward_reads_R1.fasta > t12_forward_reads_hpc_R1.fasta
seqtk hpc t12_reverse_reads_R2.fasta > t12_reverse_reads_hpc_R2.fasta

rm t12_forward_reads_R1.fasta t12_reverse_reads_R2.fasta


cat 22CMMHLT3_19570457_S63_L003_R1_001.fastq.gz 22CMMHLT3_19570459_S64_L003_R1_001.fastq.gz > t6_forward_reads_R1.fastq.gz

cat 22CMMHLT3_19570457_S63_L003_R2_001.fastq.gz 22CMMHLT3_19570459_S64_L003_R2_001.fastq.gz > t6_reverse_reads_R2.fastq.gz


seqtk seq -A t6_forward_reads_R1.fastq.gz > t6_forward_reads_R1.fasta

seqtk seq -A t6_reverse_reads_R2.fastq.gz > t6_reverse_reads_R2.fasta

seqtk hpc t6_forward_reads_R1.fasta > t6_forward_reads_hpc_R1.fasta
seqtk hpc t6_reverse_reads_R2.fasta > t6_reverse_reads_hpc_R2.fasta

rm t6_forward_reads_R1.fasta t6_reverse_reads_R2.fasta

cat 22CMMHLT3_19570451_S60_L003_R1_001.fastq.gz 22CMMHLT3_19570453_S61_L003_R1_001.fastq.gz 22CMMHLT3_19570455_S62_L003_R1_001.fastq.gz > t0_forward_reads_R1.fastq.gz

cat 22CMMHLT3_19570451_S60_L003_R2_001.fastq.gz 22CMMHLT3_19570453_S61_L003_R2_001.fastq.gz 22CMMHLT3_19570455_S62_L003_R2_001.fastq.gz > t0_reverse_reads_R2.fastq.gz


seqtk seq -A t0_forward_reads_R1.fastq.gz > t0_forward_reads_R1.fasta

seqtk seq -A t0_reverse_reads_R2.fastq.gz > t0_reverse_reads_R2.fasta

seqtk hpc t0_forward_reads_R1.fasta > t0_forward_reads_hpc_R1.fasta
seqtk hpc t0_reverse_reads_R2.fasta > t0_reverse_reads_hpc_R2.fasta

rm t0_forward_reads_R1.fasta t0_reverse_reads_R2.fasta


#Run meryl-lookup to identify strains specific reads
meryl-lookup -existence -sequence t24_forward_reads_hpc_R1.fasta  -mers /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/AB.only.meryl/ /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/TU.only.meryl/ /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/TL.only.meryl/ /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/WIK.only.meryl/ -output t24_forward_readsKemrs.txt 


meryl-lookup -existence -sequence t12_forward_reads_hpc_R1.fasta -mers /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/AB.only.meryl/ /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/TU.only.meryl/ /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/TL.only.meryl/ /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/WIK.only.meryl/ -output t12_forward_readsKemrs.txt


meryl-lookup -existence -sequence  t6_forward_reads_hpc_R1.fasta  -mers /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/AB.only.meryl/ /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/TU.only.meryl/ /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/TL.only.meryl/ /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/WIK.only.meryl/ -output t6_forward_readsKemrs.txt


meryl-lookup -existence -sequence  t0_forward_reads_hpc_R1.fasta  -mers /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/AB.only.meryl/ /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/TU.only.meryl/ /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/TL.only.meryl/ /data/okendojo/zebrafish/data/g3/merylDBs/compressedMeryls/WIK.only.meryl/ -output t0_forward_readsKemrs.txt



cat  t24_forward_readsKemrs.txt  | awk '{SUM=$4+$6+$8+$10; tag=$4":"$6":"$8":"$10; if (SUM == 0) {NAME="UNKNOWN"; color="#AAAAAA";} else if ($4/SUM > 0.9) {NAME="AB"; color="#d7191c";} else if ($6/SUM > 0.9) {NAME="TU"; color="fdae61"; } else if ($8/SUM>0.9) {NAME="TL"; color="#abdda4"; } else if ($10/SUM>0.9) { NAME="WIK"; color="#2b83ba"; } else { NAME="MIXED"; color="#FFFF00"; } print $1"\t"$2"\t"NAME"\t"color"\t"tag; }' >> t24_forward_Tags.csv

cat t12_forward_readsKemrs.txt  | awk '{SUM=$4+$6+$8+$10; tag=$4":"$6":"$8":"$10; if (SUM == 0) {NAME="UNKNOWN"; color="#AAAAAA";} else if ($4/SUM > 0.9) {NAME="AB"; color="#d7191c";} else if ($6/SUM > 0.9) {NAME="TU"; color="fdae61"; } else if ($8/SUM>0.9) {NAME="TL"; color="#abdda4"; } else if ($10/SUM>0.9) { NAME="WIK"; color="#2b83ba"; } else { NAME="MIXED"; color="#FFFF00"; } print $1"\t"$2"\t"NAME"\t"color"\t"tag; }' >> t12_forward_Tags.csv


cat t6_forward_readsKemrs.txt | awk '{SUM=$4+$6+$8+$10; tag=$4":"$6":"$8":"$10; if (SUM == 0) {NAME="UNKNOWN"; color="#AAAAAA";} else if ($4/SUM > 0.9) {NAME="AB"; color="#d7191c";} else if ($6/SUM > 0.9) {NAME="TU"; color="fdae61"; } else if ($8/SUM>0.9) {NAME="TL"; color="#abdda4"; } else if ($10/SUM>0.9) { NAME="WIK"; color="#2b83ba"; } else { NAME="MIXED"; color="#FFFF00"; } print $1"\t"$2"\t"NAME"\t"color"\t"tag; }' >> t6_forward_Tags.csv


cat t0_forward_readsKemrs.txt | awk '{SUM=$4+$6+$8+$10; tag=$4":"$6":"$8":"$10; if (SUM == 0) {NAME="UNKNOWN"; color="#AAAAAA";} else if ($4/SUM > 0.9) {NAME="AB"; color="#d7191c";} else if ($6/SUM > 0.9) {NAME="TU"; color="fdae61"; } else if ($8/SUM>0.9) {NAME="TL"; color="#abdda4"; } else if ($10/SUM>0.9) { NAME="WIK"; color="#2b83ba"; } else { NAME="MIXED"; color="#FFFF00"; } print $1"\t"$2"\t"NAME"\t"color"\t"tag; }' >> t0_forward_Tags.csv 




