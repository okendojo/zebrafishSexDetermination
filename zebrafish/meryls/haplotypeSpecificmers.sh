#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=G3_meryl

#Load module
module load meryl

cd /data/okendojo/zebrafish/data/g3/illumina/meryldbs
# AB only
meryl output AB.only.meryl difference [ difference [ difference AB.meryl TL.meryl ] TU.meryl ] WIK.meryl

# TL only
meryl output TL.only.meryl difference [ difference [ difference TL.meryl AB.meryl ] TU.meryl ] WIK.meryl

# TU only
meryl output TU.only.meryl difference [ difference [ difference TU.meryl WIK.meryl ] AB.meryl ] TL.meryl

# WIK only
meryl output WIK.only.meryl difference [ difference [ difference WIK.meryl TU.meryl ] AB.meryl ] TL.meryl

#plot shared and specific k-mers
for HAP in TU AB TL WIK ; do echo $HAP ; meryl output $HAP.shared.meryl difference $HAP.meryl $HAP.only.meryl ; $MERQURY/plot/to_hist_for_plotting.sh $HAP.shared.meryl $HAP.shared $HAP.only.meryl $HAP.only > $HAP.shared_and_only.hist ; Rscript $MERQURY/plot/plot_spectra_cn.R -f $HAP.shared_and_only.hist -o $HAP.shared_and_only -t stack -m 135 -n 35454358 ; done 
