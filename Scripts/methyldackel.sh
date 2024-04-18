#!/bin/bash

for i in $(seq 1 $1)
do
    if compgen -G "${2}_rep${i}_*bam" > /dev/null; then
        if compgen -G "${2}_rep${i}_library_sorted.bam" > /dev/null; then
            MethylDackel extract --mergeContext -l CpGIslands_mod.bed hg38.fa ${2}_rep${i}_library_sorted.bam
            MethylDackel extract --mergeContext -l CpGIslands_mod.bed --fraction hg38.fa ${2}_rep${i}_library_sorted.bam
        else
            mv ${2}_rep${i}*.bam ${2}_rep${i}_library_sorted.bam
            MethylDackel extract --mergeContext -l CpGIslands_mod.bed hg38.fa ${2}_rep${i}_library_sorted.bam
            MethylDackel extract --mergeContext -l CpGIslands_mod.bed --fraction hg38.fa ${2}_rep${i}_library_sorted.bam
        fi
    fi
    
done