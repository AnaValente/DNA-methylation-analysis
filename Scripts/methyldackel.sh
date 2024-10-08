#!/bin/bash

check_bam_file() {
    compgen -G "$1" > /dev/null
}

process_bam_file() {
    local bam_file=$1
    local prefix=$2
    local library_sorted_bam="${prefix}_library_sorted.bam"
    local methyl
    local condition

    # Check if library_sorted.bam exists, if not, rename existing bam file
    if ! check_bam_file "${library_sorted_bam}"; then
        mv ${prefix}*.bam ${library_sorted_bam}
    fi

    # Run MethylDackel mbias and extract condition
    if [ "$3" == "false" ]; then
        methyl="$(MethylDackel mbias genome.fa ${library_sorted_bam} ${prefix} 2>&1)"
    else
        methyl="$(MethylDackel mbias -l $3 genome.fa ${library_sorted_bam} ${prefix} 2>&1)"
    fi

    condition="$(echo ${methyl} | sed -e 's/^.*://g')"

    # Run MethylDackel extract
    if [ "$3" == "false" ]; then
        MethylDackel extract --mergeContext genome.fa ${condition} ${library_sorted_bam}
    else
        MethylDackel extract --mergeContext -l $3 genome.fa ${condition} ${library_sorted_bam}
    fi
}

for i in $(seq 1 $1); do
    prefix="${2}_rep${i}"

    if check_bam_file "${prefix}*bam"; then
        process_bam_file "${prefix}" "${prefix}" "$3"
    fi
done