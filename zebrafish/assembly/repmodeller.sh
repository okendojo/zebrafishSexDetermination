#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=repmodeller


#Load the modules
module load repeatmasker
module load repeatmodeler
module load trf

cd /data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/centromere

# ==== CONFIG ====
GENOME="/data/okendojo/zebrafish/data/fish6/asm/t2t/polishing/asm_polished/NHGRI_Fish6_cons.fasta"        # Input genome FASTA
DBNAME="zebrafish_db"             # RepeatModeler DB name
THREADS=16                         # Number of threads
OUTDIR="repeat_output"            # Output directory

# ==== SETUP ====
mkdir -p $OUTDIR
cd $OUTDIR

echo "Step 1: Building RepeatModeler database..."
BuildDatabase -name $DBNAME $GENOME

echo "Step 2: Running RepeatModeler (this may take hours)..."
RepeatModeler -database $DBNAME -pa $THREADS

echo "Step 3: Running RepeatMasker with custom library..."
RepeatMasker -pa $THREADS -lib ${DBNAME}-families.fa -gff -dir masked $GENOME

echo "Step 4: Running Tandem Repeat Finder..."
trf $GENOME 2 7 7 80 10 50 2000 -d -h

echo "Step 5: Searching for candidate centromeric repeats..."
grep -iE 'sat|cen|alpha|hor' ${DBNAME}-families.fa > candidate_centromeric_repeats.fa
grep -iE 'sat|cen|alpha|hor' masked/${GENOME}.out > candidate_repeats_masked.out

echo "✅ Done. Outputs:"
echo "- Custom repeat library: ${DBNAME}-families.fa"
echo "- Masked genome: masked/${GENOME}.masked"
echo "- TRF output: ${GENOME}.2.7.7.80.10.50.2000.dat"
echo "- Candidate repeats: candidate_centromeric_repeats.fa"

