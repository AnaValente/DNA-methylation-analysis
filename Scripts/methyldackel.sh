#!/bin/bash

for i in $(seq 1 $1)
do
    if [ -f ${2}_rep2_library_sorted.bam ]; then
        if [ -f ${2}_rep${i}_library_sorted.bam ]; then
            MethylDackel extract --mergeContext -l CpGIslands_mod.bed hg38.fa ${2}_rep${i}_library_sorted.bam
            MethylDackel extract --mergeContext -l CpGIslands_mod.bed --fraction hg38.fa ${2}_rep${i}_library_sorted.bam
        else
            mv ${2}_rep${i}*.bam ${2}_rep${i}_library_sorted.bam
            MethylDackel extract --mergeContext -l CpGIslands_mod.bed hg38.fa ${2}_rep${i}_library_sorted.bam
            MethylDackel extract --mergeContext -l CpGIslands_mod.bed --fraction hg38.fa ${2}_rep${i}_library_sorted.bam
        fi
    else
        if [ -f ${2}_rep${i}_library_sorted.bam ]; then
            MethylDackel extract --mergeContext -l CpGIslands_mod.bed hg38.fa ${2}_rep1_library_sorted.bam
            MethylDackel extract --mergeContext -l CpGIslands_mod.bed --fraction hg38.fa ${2}_rep1_library_sorted.bam
        else
            mv ${2}_rep1*.bam ${2}_rep1_library_sorted.bam
            MethylDackel extract --mergeContext -l CpGIslands_mod.bed hg38.fa ${2}_rep1_library_sorted.bam
            MethylDackel extract --mergeContext -l CpGIslands_mod.bed --fraction hg38.fa ${2}_rep1_library_sorted.bam
        fi
    
    fi
    
done