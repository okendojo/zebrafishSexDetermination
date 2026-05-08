#! /bin/sh

###############
##HEADER FOR BIOWULF SBATCH SUBMIT
##############
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --partition=norm
#SBATCH --mail-type=FAIL,END
#SBATCH --mem=100G
#SBATCH --gres=lscratch:400
#SBATCH --time=240:00:00
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=BRAKER
#######
## LOAD THE MODULES NEEDED TO RUN THE WORKFLOW
#######
module load braker/3
module add diamond
module add blast
module load python
module load perl

#export GENEMARK_PATH=/home/okendojo/software/geneMark/gmes_linux_64
#export AUGUSTUS_BIN_PATH=/home/okendojo/software/config

#export AUGUSTUS_BIN_PATH=/data/okendojo/conda/envs/BRAKER/bin/


########
## DEFINE FILES AND FILE LOCATIONS
########

## number of cores to use; set THREADS to a multiple of 4; should match SBATCH
## max memory to use; should match SBATCH
THREADS=32
MEM=100G



##HOME_DIR is where you will direct outputs

HOME_DIR="/data/okendojo/zebrafish/data/fish11/annotation"
##GENOME is the source of the genome FASTA

GENOME="/data/Zebrafish_T2T/fish11/resolved_assembly/fish11_merfin_polished.fasta"
PROTSEQ="/data/okendojo/zebrafish/data/fish11/annotation/uniprotkb_zebrafish_2023_10_05.fasta"
PROHINT="/data/okendojo/zebrafish/data/fish11/annotation/protHint/prothint_augustus.gff"
##SPECIES_NAME should be the name you want the project files to be under

SPECIES_NAME="zebrafish"
	
#############
## script commands below
#############

cd ${HOME_DIR}



## BRAKER
## Running the first round of braker with 


echo "----- beginning braker 2 at" `date`

#braker.pl --threads=$THREADS --verbosity=3  --species=$SPECIES_NAME --AUGUSTUS_BIN_PATH="/data/okendojo/conda/envs/BRAKER/bin" --AUGUSTUS_CONFIG_PATH="/data/okendojo/conda/envs/BRAKER/config" --AUGUSTUS_SCRIPTS_PATH="/usr/local/apps/augustus/3.4.0/scripts" --genome=$GENOME --esmode --useexisting



#Copy of the origian run
braker.pl --threads=$THREADS --verbosity=3  --species=$SPECIES_NAME --AUGUSTUS_BIN_PATH="/data/okendojo/conda/envs/BRAKER/bin" --AUGUSTUS_CONFIG_PATH="/data/okendojo/conda/envs/BRAKER/config" --AUGUSTUS_SCRIPTS_PATH="/usr/local/apps/augustus/3.4.0/scripts" --genome=$GENOME --BAMTOOLS_PATH=/data/okendojo/conda/envs/BRAKER/bin/ --prot_seq=$PROTSEQ --useexisting --rnaseq_sets_ids="ERR10034048,ERR10034049,ERR10034050,ERR10034051,ERR10034052,ERR10034068,ERR10034069,ERR10034070,ERR10034071,ERR10034072,ERR10034088,ERR10034089,ERR10034090"

echo "----- finished braker_2 at" `date`

