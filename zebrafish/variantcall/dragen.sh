#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=64g
#SBATCH --ntasks-per-core=1
#SBATCH --time=96:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=javan.okendo@nih.gov

#load the module
module load bwa
module load samtools
module load dragmap/1.2.1
module load GATK
module load picard

cd /data/okendojo/paradisfishProject/snv/dragmapping

#map the reads
#dragen-os -r /data/okendojo/paradisfishProject/snv/dragenRef -1 ../ERR3332352_1.fastq.gz -2 ../ERR3332352_2.fastq.gz > macOpe2.sam

#get str table
#gatk ComposeSTRTableFile -R ../assembly/macOpe2Assembly.fasta -O str_table.tsv

#convert sam to bam
#samtools view -S -b macOpe2.sam -o macOpe2.bam
#samtools sort macOpe2.bam -o macOpe2.sorted.bam
#samtools index macOpe2.sorted.bam

#Add read groups else haplotype caller will not run
#java -jar $PICARDJARPATH/picard.jar AddOrReplaceReadGroups -I macOpe2.bam  -O macOpe2_RG.bam -RGID 4 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM 20 -SO coordinate --CREATE_INDEX true

#calibrate model
#gatk CalibrateDragstrModel \
 #   -R ../assembly/macOpe2Assembly.fasta \
  #  -I macOpe2.sorted.bam \
   # -str str_table.tsv \
    #-O dragstr_model.txt

#call variants
gatk HaplotypeCaller \
    -R ../assembly/macOpe2Assembly.fasta \
    -I macOpe2.sorted.bam \
    -O sv_output_file.vcf \
    --dragen-mode true \
    --add-output-vcf-command-line false \
    --dragstr-params-path dragstr_model.txt

#filter variants
#gatk VariantFiltration \
 #     -V output_file.vcf \
  #    --filter-expression "QUAL < 10.4139" \
   #   --filter-name "DRAGENHardQUAL" \
    #  -O output_filtered.vcf





