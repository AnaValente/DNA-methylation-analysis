#!/bin/bash

for i in $(seq 1 $1)
do
    if compgen -G "${2}_rep${i}*bam" > /dev/null; then
        if compgen -G "${2}_rep${i}_library_sorted.bam" > /dev/null; then
            methyl="$(MethylDackel mbias -l CpGIslands_mod.bed genome.fa ${2}_rep${i}_library_sorted.bam ${2}_rep${i} 2>&1)"
            condition="$(echo ${methyl} | sed -e 's/^.*://g')"
            MethylDackel extract --mergeContext -l CpGIslands_mod.bed genome.fa ${condition} ${2}_rep${i}_library_sorted.bam
        else
            mv ${2}_rep${i}*.bam ${2}_rep${i}_library_sorted.bam
            methyl="$(MethylDackel mbias -l CpGIslands_mod.bed genome.fa ${2}_rep${i}_library_sorted.bam ${2}_rep${i} 2>&1)"
            condition="$(echo ${methyl} | sed -e 's/^.*://g')"
            MethylDackel extract --mergeContext -l CpGIslands_mod.bed genome.fa ${condition} ${2}_rep${i}_library_sorted.bam
        fi
    fi
    
done