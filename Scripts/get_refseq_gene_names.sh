#!/usr/bin/env bash 

# search for closest coordinates between RefSeq and DMRs and store in bedgraph with gene names
sort-bed $1 | bedtools closest -D ref -a - -b  $2 | cut -f1,2,3,4,8,9 | awk '!x[$0]++' - > $3_refseq.bedgraph