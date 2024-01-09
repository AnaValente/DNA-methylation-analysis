#!/usr/bin/env bash 

# download RefSeq genes
curl -s "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/ncbiRefSeq.txt.gz" | gunzip -c | cut -f3,5,6,13 | sort-bed - > genes.bed

# search for overlaping coordinates between RefSeq and DMRs and store in bedgraph with gene names
sort-bed $1* | bedtools closest -D ref -a - -b genes.bed | cut -f1,2,3,4,8,9 | awk '!x[$0]++' - > $2_refseq.bedgraph


