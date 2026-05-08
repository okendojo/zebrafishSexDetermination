#!/bin/bash
set -e -o pipefail
module load python/3.8
module load singularity
module load pytorch
module load pigz
export PATH=$PATH:/home/korens/.local/bin/

#  Discover the job ID to run, from either a grid environment variable and a
#  command line offset, or directly from the command line.
#
if [ x$SLURM_ARRAY_TASK_ID = x -o x$SLURM_ARRAY_TASK_ID = xundefined -o x$SLURM_ARRAY_TASK_ID = x0 ]; then
  baseid=$1
  offset=0
  cores=`grep -c ^processor /proc/cpuinfo`
else
  baseid=$SLURM_ARRAY_TASK_ID
  offset=$1
  cores=$SLURM_CPUS_PER_TASK
fi
if [ x$offset = x ]; then
  offset=0
fi
if [ x$baseid = x ]; then
  echo Error: I need SLURM_ARRAY_TASK_ID set, or a job index on the command line.
  exit 0
fi
jobid=`expr -- $baseid + $offset |awk '{printf("%03d", $1)}'`
if [ x$SLURM_ARRAY_TASK_ID = x ]; then
  echo Running job $jobid based on command line options.
else
  echo Running job $jobid based on SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID and offset=$offset.
fi

MAX_ID=`find -L -maxdepth 1 -type d |grep -v "^\.$"|grep -v github | wc -l`
if [ x$MAX_ID = x ]; then
   MAX_ID=0
fi

if [ $jobid -gt $MAX_ID ]; then
  echo Error: Only $MAX_ID partitions, you asked for $jobid.
  exit 0
fi

sample=`find -L -maxdepth 1 -type d |grep -v "^\.$" |grep -v github | sort | head -n $jobid |tail -n 1`
echo "Running $jobid with $sample"

if [ ! -e $sample/herro/corrected.fasta.gz ]; then
   mkdir -p $sample/herro
   cd $sample/herro
   ONTinput=`ls ../*.bam |grep -i 5mc`
   isR9=`samtools view -H $ONTinput |grep -c r9 || :`
   isR10=`samtools view -H $ONTinput |grep -c r10 || :`
   if [ $isR9 -ge 1 ] && [ $isR10 -ge 1 ]; then
      echo "Mixed R9 and R10 sequence type"
	  exit 1
   elif [ $isR9 -ge 1 ]; then
      model="/data/korens/devel/herro/model/herro_R9-v1"
   elif [ $isR10 -ge 1 ]; then
      model="/data/korens/devel/herro/model/herro-v1"
   else
      echo "Unknown sequence type, exiting"
      exit 1
   fi
   echo "Starting with model $model from $ONTinput"
exit
   F=`ls ../*readlen |grep -i 5mc`
   MIN_LEN=`cat $F | sort -srnk2,2 | awk -v G=3100000000 'BEGIN { L=0; } {S+=$2; if (S/G >= 40 && L == 0) L=$2; } END {print L}'`
   echo "Subsampling to reads longer than $MIN_LEN (from $F) for 40x and correcting with model $model"

   if [ ! -e input.fastq ]; then
      > input.fastq
      for j in `ls $ONTinput` ; do
	     echo "Adding reads from $j"
         samtools fastq $j | seqtk seq -L $MIN_LEN - >> input.fastq || exit 1
      done
      cat $F |awk -v M=$MIN_LEN '{if ($2 >= M) print $3}' |grep -v sequence-name > input.ids
   fi

   module load dorado
   module load seqkit
   module load porechop

   res="cpu"
   icores=$cores
   batches=8

   if [ "x$CUDA_VISIBLE_DEVICES" = "x" ]; then
	  res="cpu"
	  icores=$cores
	  batches=$cores
   else
      # check gpu mem and set batches accordingly
      mem=`nvidia-smi -L|awk '{print $4}'|sort |uniq |awk -F "-" '{print $NF}' |sed s/GB//g`

      if [ $mem -le 24 ]; then
         batches=2
      elif [ $mem -gt 64 ]; then
         batches=40
      fi
      icores=`echo "$cores" |awk -v S=$SLURM_GPUS_ON_NODE '{print $1/S}'`
      res="cuda:$CUDA_VISIBLE_DEVICES"
   fi
   echo "Running herro on device $res w/$SLURM_GPUS_ON_NODE and so $icores out of $SLURM_CPUS_PER_TASK and GPU devices $CUDA_VISIBLE_DEVICES"
   dorado correct -t $cores -x $res --infer-threads $icores -b $batches -i 4G -l input.ids -m $model input.fastq 2> dorado.err > corrected.fasta && pigz corrected.fasta && rm ./input.*
fi
