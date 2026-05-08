#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=230g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=hapmers

cd /data/okendojo/zebrafish/data/g3/merylDBs/f1gen/hapmers 

$MERQURY/trio/hapmers.sh ../mat/maternal.k30.meryl ../pat/paternal.k30.meryl ../child/child.k30.meryl


