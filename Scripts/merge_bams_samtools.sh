#!/bin/bash

# merge bams of different runs
samtools merge $1_library.bam $1*.bam && samtools sort $1_library.bam -o results/$1_library_sorted.bam -O BAM &&  samtools index results/$1_library_sorted.bam

rm $1_library.bam 
