#!/bin/bash
#This work was adopted from: https://bioinformaticsworkbook.org/dataWrangling/R/visualize-gaps-in-genomes.html#gsc.tab=0
genome="$1"
out=$(basename ${genome%.*})
module load bioawk
tr "ATGCatgc" "xxxxxxxx" < $genome |\
bioawk -c fastx '{print $seq}' | \
sed -e 's/xN/\nN/g' -e 's/Nx/N\n/g' | \
sed 's/x//g' | \
awk '{print length($1)}' |\
grep -vw "^0$" > ${out}_Nsizes.txt
