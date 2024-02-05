#!/bin/bash

samples_rep=""
control_rep=""

if  [[ ! -f ${1}_rep2_library_sorted_CpG.meth.bedGraph ]] && [[ -f ${2}_rep2_library_sorted_CpG.meth.bedGraph ]]; then
    metilene_input.pl --in1 ${2}_rep1_library_sorted_CpG.meth.bedGraph,${2}_rep2_library_sorted_CpG.meth.bedGraph --in2 ${1}_rep1_library_sorted_CpG.meth.bedGraph
    metilene -a g1 -b g2 metilene_g1_g2.input | sort -V -k1,1 -k2,2n > metilene_${2}_${1}.bedgraph
    metilene_output.pl -q metilene_${2}_${1}.bedgraph -o ./metilene_${2}_${1}

else

    for i in $(seq 1 ${3})
    do
        samples_rep=${samples_rep}",${2}_rep${i}_library_sorted_CpG.meth.bedGraph"
        control_rep=${control_rep}",${1}_rep${i}_library_sorted_CpG.meth.bedGraph"
    done

    metilene_input.pl --in1 ${samples_rep/,/} --in2 ${control_rep/,/}
    metilene -a g1 -b g2 metilene_g1_g2.input | sort -V -k1,1 -k2,2n > metilene_${2}_${1}.bedgraph
    metilene_output.pl -q metilene_${2}_${1}.bedgraph -o ./metilene_${2}_${1}

fi
