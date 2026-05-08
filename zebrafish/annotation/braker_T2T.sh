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
#SBATCH --job-name=BRAKER_T2T
#######
## LOAD THE MODULES NEEDED TO RUN THE WORKFLOW
#######
module load braker/3
module add diamond
module add java
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

##SPECIES_NAME should be the name you want the project files to be under

SPECIES_NAME="zebrafish"
	
#############
## script commands below
#############

cd ${HOME_DIR}



## BRAKER
## Running the first round of braker with 


echo "----- beginning braker3 at" `date`


#Copy of the origian run
braker.pl --threads=$THREADS --verbosity=3  --species=$SPECIES_NAME  --workingdir="/data/okendojo/zebrafish/data/fish11/annotation/T2T_annotation" --AUGUSTUS_BIN_PATH="/data/okendojo/conda/envs/BRAKER/bin" --AUGUSTUS_CONFIG_PATH="/data/okendojo/conda/envs/BRAKER/config" --AUGUSTUS_SCRIPTS_PATH="/usr/local/apps/augustus/3.4.0/scripts" --genome=$GENOME --BAMTOOLS_PATH=/data/okendojo/conda/envs/BRAKER/bin/ --prot_seq=$PROTSEQ  --useexisting --rnaseq_sets_ids="SRR24844105,SRR24844106,SRR24844107,SRR24844108,SRR24844109,SRR24844110,SRR24844111,SRR24844112,SRR24844113,SRR24844114,SRR24042928,SRR24042929,SRR24042930,SRR24042931,SRR24042932,SRR24042933,SRR24042934,SRR24042935,SRR24042936,SRR24042937,SRR24057966,SRR24057967,SRR24057968,SRR24057969,SRR24057970,SRR24057971,SRR24057972,SRR24057973,SRR25232723,SRR25232724,SRR25232725,SRR25232726,SRR25232727,SRR25232728,SRR25232729,SRR25232730,SRR25232731,SRR25232732,ERR10297495,ERR10306120,ERR10306121,ERR10306122,ERR10306123,ERR10306124,ERR10306125,ERR10306126,ERR10306127,ERR10306128,ERR11282336,ERR11282337,ERR11282339,ERR11282340,ERR11282341,ERR11282342,ERR11282343,ERR11282344,ERR11282345,ERR11282346,SRR26321412,SRR26321413,SRR26321414,SRR26321415,SRR26321416,SRR26321417,SRR26321418,SRR26321419,SRR26321420,SRR26321421"

echo "----- finished braker3 at" `date`

