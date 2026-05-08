#! /bin/sh

###############
##HEADER FOR BIOWULF SBATCH SUBMIT
##############

#SBATCH --job-name=T2T_Annotation
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --partition=norm
#SBATCH --mail-type=FAIL,END
#SBATCH --mem=100G
#SBATCH --gres=lscratch:400
#SBATCH --time=240:00:00
#SBATCH --mail-user=javan.okendo@nih.gov

#######
## ACTIVATE mamba ENVIRONMENT
#######

#move to the dir containing bam files

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

HOME_DIR="/data/okendojo/zebrafish/data/fish11/annotation/T2T_annotation/rnaseq/hisat2"

##GENOME is the source of the genome FASTA

GENOME="/data/Zebrafish_T2T/fish11/resolved_assembly/fish11_merfin_polished.fasta"

## Transcripts aligned to genome in BAM format

TRANSCRIPT_BAM="ERR10297495.bam,  ERR10306125.bam  ,ERR11282339.bam  ,ERR11282345.bam  ,SRR24042932.bam , SRR24057966.bam  ,SRR24057972.bam  ,SRR24844109.bam , SRR25232723.bam , SRR25232729.bam , SRR26321414.bam,  SRR26321420.bam,
ERR10306120.bam  ,ERR10306126.bam , ERR11282340.bam , ERR11282346.bam , SRR24042933.bam , SRR24057967.bam,  SRR24057973.bam , SRR24844110.bam , SRR25232724.bam,  SRR25232730.bam,  SRR26321415.bam , SRR26321421.bam,
ERR10306121.bam  ,ERR10306127.bam , ERR11282341.bam , SRR24042928.bam , SRR24042934.bam , SRR24057968.bam , SRR24844105.bam , SRR24844111.bam  ,SRR25232725.bam , SRR25232731.bam , SRR26321416.bam,
ERR10306122.bam , ERR10306128.bam , ERR11282342.bam  ,SRR24042929.bam  ,SRR24042935.bam  ,SRR24057969.bam , SRR24844106.bam , SRR24844112.bam  ,SRR25232726.bam , SRR25232732.bam  ,SRR26321417.bam,
ERR10306123.bam  ,ERR11282336.bam  ,ERR11282343.bam  ,SRR24042930.bam  ,SRR24042936.bam  ,SRR24057970.bam  ,SRR24844107.bam , SRR24844113.bam  ,SRR25232727.bam , SRR26321412.bam  ,SRR26321418.bam,
ERR10306124.bam  ,ERR11282337.bam  ,ERR11282344.bam  ,SRR24042931.bam  ,SRR24042937.bam  ,SRR24057971.bam  ,SRR24844108.bam  ,SRR24844114.bam  ,SRR25232728.bam  ,SRR26321413.bam  ,SRR26321419.bam"

##SPECIES_NAME should be the name you want the project files to be under

SPECIES_NAME="zebrafish"
	
#############
## script commands below
#############

cd ${HOME_DIR}

## BRAKER
## Running the first round of braker with 


echo "----- beginning braker 2 at" `date`

braker.pl --threads=$THREADS --AUGUSTUS_BIN_PATH="/data/okendojo/conda/envs/BRAKER/bin" --AUGUSTUS_CONFIG_PATH="/data/okendojo/conda/envs/BRAKER/config" --AUGUSTUS_SCRIPTS_PATH="/usr/local/apps/augustus/3.4.0/scripts" --verbosity=4 --makehub --email "javanokendo@gmail.com" --species=$SPECIES_NAME --genome=$GENOME --bam=$TRANSCRIPT_BAM --BAMTOOLS_PATH="/data/okendojo/conda/envs/BRAKER/bin/" --softmasking --useexisting 

echo "----- finished braker_2 at" `date`

