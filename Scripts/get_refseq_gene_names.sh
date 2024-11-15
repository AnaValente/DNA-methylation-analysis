#!/usr/bin/env bash 

# Search for closest coordinates between RefSeq file and DMR file and store in bedgraph with gene names
sort-bed $1 | bedtools closest -D ref -a - -b  $2 | cut -f1,2,3,4,8,9 | awk '!x[$0]++' - > $3_refseq.bedgraph