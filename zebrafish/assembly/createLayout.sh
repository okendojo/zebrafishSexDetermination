#!/bin/sh
set -e

 /usr/local/apps/verkko/2.1/lib/verkko/scripts/get_layout_from_mbg.py \
  /data/okendojo/zebrafish/data/fish6/asm/t2t_chrs/6-layoutContigs/combined-nodemap.txt \
  /data/okendojo/zebrafish/data/fish6/asm/t2t_chrs/6-layoutContigs/combined-edges.gfa \
  /data/okendojo/zebrafish/data/fish6/asm/t2t_chrs/6-layoutContigs/combined-alignments.gaf \
  /data/okendojo/zebrafish/data/fish6/asm/consensusfile/tmp.txt \
  /data/okendojo/zebrafish/data/fish6/asm/t2t_chrs/6-layoutContigs/nodelens.txt \
  /data/okendojo/zebrafish/data/fish6/asm/t2t_chrs/6-layoutContigs/unitig-popped.layout \
  /data/okendojo/zebrafish/data/fish6/asm/t2t_chrs/6-layoutContigs/unitig-popped.layout.scfmap

 /usr/local/apps/verkko/2.1/lib/verkko/scripts/check_layout_gaps.py < /data/okendojo/zebrafish/data/fish6/asm/t2t_chrs/6-layoutContigs/unitig-popped.layout > /data/okendojo/zebrafish/data/fish6/asm/t2t_chrs/6-layoutContigs/gaps.txt
