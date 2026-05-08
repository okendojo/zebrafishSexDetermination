#!/bin/bash
#SBATCH --partition=norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=180g
#SBATCH --ntasks-per-core=1
#SBATCH --time=240:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=javan.okendo@nih.gov
#SBATCH --job-name=MGG_fish11


#Load the modules
module load mbg

cd /data/okendojo/zebrafish/data/fish11/asm_vk210x_v1_4

MBG -i hifi-corrected.fasta.gz -o hifi_graph.gfa -k 1501 -w 1450 -a 1 -u 3 -t 24 -r 15000 -R 4000 

awk 'BEGIN \
     { \
        FS="[ \t]+"; OFS="\t"; \
        print "node", "length", "coverage"; \
     } \
     $1 == "S" \
     { \
        if ($6 != "") {
          $4 = $6;
        }
        print $2, length($3), substr($4, 6); \
     }' hifi_graph.gfa > hifi_graph.gfa_nodecov.csv
