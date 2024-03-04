#!/usr/bin/env nextflow

//params
params.files = false
params.samples = false
params.concat = false
params.replicates = false
params.cell_tpm = false
params.cutoff_regions = 75
params.cutoff_heatmap = 100
params.genome = false

c_green = "\033[0;38;5;28m"
c_yellow = "\033[0;38;5;226m"
c_dark_blue = "\033[38;5;75m"
c_orange = "\033[38;5;202m"
c_blue = "\033[0;38;5;33m"
c_reset = "\033[0m"

log.info"""
       ${c_dark_blue}. ,-"-.   ,-"-. ,-"-.   ,-"-. ,-"-.   ,-"-. ,-"-.    ,-"-. ,-"-.   ,${c_reset}
        ${c_dark_blue}/ ${c_green}| ${c_yellow}| ${c_dark_blue}\\ / ${c_orange}| ${c_blue}| ${c_dark_blue}/ ${c_blue}| ${c_green}| ${c_dark_blue}\\ / ${c_yellow}| ${c_blue}| ${c_dark_blue}/ ${c_yellow}| ${c_orange}| ${c_dark_blue}\\ / ${c_green}| ${c_blue}| ${c_dark_blue}/ ${c_orange}| ${c_blue}| ${c_dark_blue}\\  / ${c_yellow}| ${c_green}| ${c_dark_blue}/ ${c_orange}| ${c_blue}| ${c_dark_blue}\\ /${c_reset}
        D N A   M E T H Y L A T I O N   A N A L Y S I S   P I P E L I N E
       ${c_dark_blue}/ \\${c_blue}| ${c_orange}| ${c_dark_blue}|\\| ${c_yellow}| ${c_green}|${c_dark_blue}/ \\${c_green}| ${c_blue}| ${c_dark_blue}|\\| ${c_orange}| ${c_green}|${c_dark_blue}/ \\${c_orange}| ${c_yellow}| ${c_dark_blue}|\\| ${c_blue}| ${c_green}|${c_dark_blue}/ \\${c_yellow}| ${c_green}| ${c_dark_blue}|\\|  ${c_orange}| ${c_blue}|${c_dark_blue}/ \\${c_yellow}| ${c_green}| ${c_dark_blue}|\\${c_reset}
          ${c_dark_blue}`-!-' `-!-'   `-!-' `-!-'   `-!-' `-!-'   `-!-'  `-!-'   `-!-' `-${c_reset}
        """.stripIndent()

process REFERENCE_GENOME_HG38 {
    // Download fasta file of the reference genome 
    
    publishDir "./Pipeline_results/Methyldackel", mode: "copy"
    
    input:
    path genome
    
    output:
    path "*"

    """
    gunzip -c ${genome} > hg38.fa
    """
}

process CONCATENATE_BAMS {
    // Concatenate bam files for each sample

    publishDir "./Pipeline_results/Concatenated_files", mode: "copy"
    
    input:
    path folder
    val samples
    val n_replicates

    output:
    path "${samples}*_library_sorted.b*"

    shell:
    '''
    for i in $(seq 1 !{n_replicates})
    do
        if compgen -G "!{samples}_rep${n}*.bam" > /dev/null; then
            samtools merge !{samples}_rep${n}_library.bam !{samples}_rep${n}*.bam && samtools sort !{samples}_rep${n}_library.bam -o !{samples}_rep${n}_library_sorted.bam -O BAM && samtools index !{samples}_rep${n}_library_sorted.bam
        fi
    done
   '''
}

process METHYLDACKEL {
    // Methylation calling

    publishDir './Pipeline_results/Methyldackel', mode: 'copy'
    
    input:
    path folder
    path bam_rep1
    path hg38
    val samples
    val n_replicates

    output:
    path "*.bedGraph"

    """
    ./methyldackel.sh ${n_replicates} ${samples}
    """
}

process CORRELATION_AND_CLUSTERING {
    // correlation and clustering analysis
	
    input:
    path folder
    path bedgraphs
    val samples
    val number_samples
    val control
    val cutoff_regions
    val cutoff_heatmap
    val n_replicates

    output:
    path "*"

    publishDir './Pipeline_results', pattern:"*png" , mode: 'copy'
    publishDir './Pipeline_results', pattern:"*pdf" , mode: 'copy'

    """
	python3 methyldackel_cluster_table.py ${bedgraphs} ${samples} ${number_samples} ${control} ${cutoff_regions} ${cutoff_heatmap} ${n_replicates} && ./correlation_analysis.R correlation.csv ${number_samples} ${n_replicates} ${control} && ./heatmap.R diff_methylation_filtered.csv 
    """
}

process METILENE { 
    // Calculation of the differentially methylated regions

    publishDir './Pipeline_results/Metilene', mode: 'copy'
    
	
    input:
    path folder
    path bedgraphs
    val control
    val samples
    val n_replicates

    output:
    path "*_qval.0.05.bedgraph"

    shell:
    '''
    samples_rep=""
    control_rep=""

    if  [[ ! -f !{control}_rep2_library_sorted_CpG.meth.bedGraph ]] && [[ -f !{samples}_rep2_library_sorted_CpG.meth.bedGraph ]]; then
        metilene_input.pl --in1 !{samples}_rep1_library_sorted_CpG.meth.bedGraph,!{samples}_rep2_library_sorted_CpG.meth.bedGraph --in2 !{control}_rep1_library_sorted_CpG.meth.bedGraph
        metilene -a g1 -b g2 metilene_g1_g2.input | sort -V -k1,1 -k2,2n > metilene_!{samples}_!{control}.bedgraph
        metilene_output.pl -q metilene_!{samples}_!{control}.bedgraph -o ./metilene_!{samples}_!{control}

    else

        for i in $(seq 1 !{n_replicates})
        do
            samples_rep=${samples_rep}",!{samples}_rep${i}_library_sorted_CpG.meth.bedGraph"
            control_rep=${control_rep}",!{control}_rep${i}_library_sorted_CpG.meth.bedGraph"
        done

        metilene_input.pl --in1 ${samples_rep/,/} --in2 ${control_rep/,/}
        metilene -a g1 -b g2 metilene_g1_g2.input | sort -V -k1,1 -k2,2n > metilene_!{samples}_!{control}.bedgraph
        metilene_output.pl -q metilene_!{samples}_!{control}.bedgraph -o ./metilene_!{samples}_!{control}

    fi
    '''
}

process GENOMIC_REGIONS {
    // Genomic regions retrieval

    publishDir './Pipeline_results', mode: 'copy'
    
	
    input:
    path folder
    path methyldackel
    path metilene

    output:
    path "*"

    """
    ./genomic_regions.R
    """
}

process GENE_NAMES {
    // Gene name retrieval

    publishDir './Pipeline_results/Gene_names', mode: 'copy'
    
	
    input:
    path folder
    path bedgraphs
    val control
    val samples
    val cell_tpm

    output:
    path "*_merged_genes*.*"

    """
	./get_refseq_gene_names.sh metilene_${samples}_${control}_qval.0.05.bedgraph ${samples}_${control} && python3 expressed_cell_genes.py ${control} ${samples} ${cell_tpm}
    """
}

process VENN_DIAGRAM {
    // Venn Diagram

    publishDir './Pipeline_results', mode: 'copy'
    
	
    input:
    path folder
    path genes


    output:
    path "*"

    """
	./venn_diagram.R
    """
}

process FINAL_REPORT {
    
    publishDir './Pipeline_results', mode: 'copy'

    input:
    path bams
    path methyldackel
    path metilene
    path correlation
    val cutoff_heatmap
    val samples

    output:
    path "Final_report.txt"

    shell:
    '''
    if [ -f *.bam ]; then
        echo -e "Concatenated BAMS:\\n=================" >> Final_report.txt
        for bam in *.bam; do
            samtools flagstat ${bam}
        done
    fi

    echo -e "Methyldackel:\\n=============" >> Final_report.txt
    wc -l *meth.bedGraph | sed 's/_library.*//g' | awk '{print $2 ":", $1-1, "CpG sites\\n"}' >> Final_report.txt
    echo -e "\\nCorrelation and clustering:\\n===========================" >> Final_report.txt
    echo -e "$(wc -l correlation.csv | awk '{print $1-1}')" "CpG sites were used for the correlation matrix and PCA\\n" >> Final_report.txt
    echo -e "$(wc -l diff_methylation_filtered.csv | awk '{print $1-1}')" "CpG sites were used for the hierarchical clustering heatmap with a difference of !{cutoff_heatmap}% between the controls and samples\\n" >> Final_report.txt
    echo -e "\\nMetilene:\\n=========" >> Final_report.txt
    
    for sample in !{samples}; do
        echo -ne "$(wc -l metilene_${sample}*_qval.0.05.bedgraph | sed 's/metilene_//' | sed 's/_.*//g' | awk '{print $2 ":", $1, "differentially methylated regions,"}') where $(awk '$4 > 0 {print ;}' metilene_${sample}*_qval.0.05.bedgraph  | wc -l) are hypermetilated" >> Final_report.txt
        echo -ne " and" "$(awk '$4 < 0 {print ;}' metilene_${sample}*_qval.0.05.bedgraph  | wc -l)" "are hypometilated\\n" >> Final_report.txt
    done
    '''
}

workflow {
    
    if (!params.files) {
        error("No input folder provided --files")
    }

    if (!params.samples) {
        error("No sample names provided --samples")
    }

    if (!params.genome) {
        error("No reference genome hg38 file provided --genome")
    }

    if (!params.replicates) {
        error("No replicate number provided --replicates")
    }
    
    if (params.cutoff_regions < 0 || params.cutoff_regions > 100 || params.cutoff_heatmap < 0 || params.cutoff_heatmap > 100) {
        error("Please provide a cutoff number between 0 and 100")
    }

    files = Channel.fromPath(params.files).collect()
    list_samples = params.samples.split(',') as List
    samples_chn = Channel.from(list_samples)
    control = list_samples.remove(0)
    samples_chn2 = Channel.from(list_samples)
    samples_names = list_samples.toString().replace("[", "").replace("]", "").replace(",", "")

    REFERENCE_GENOME_HG38(Channel.fromPath(params.genome))
    
    if (params.concat) {
        CONCATENATE_BAMS(files, samples_chn, params.replicates)
        METHYLDACKEL(files, CONCATENATE_BAMS.out.collect(), REFERENCE_GENOME_HG38.out.collect(), samples_chn,params.replicates)
    }
    else {
        METHYLDACKEL(files, [], REFERENCE_GENOME_HG38.out.collect(), samples_chn,params.replicates)
    }

    CORRELATION_AND_CLUSTERING(files, METHYLDACKEL.out.collect(), samples_names, list_samples.size(), control, params.cutoff_regions, params.cutoff_heatmap, params.replicates)
    METILENE(files, METHYLDACKEL.out.collect(), control, samples_chn2, params.replicates)
    GENOMIC_REGIONS(files, CORRELATION_AND_CLUSTERING.out.collect(), METILENE.out.collect())
    GENE_NAMES(files, METILENE.out.collect(), control, samples_chn2, params.cell_tpm)

    if (list_samples.size() > 1 & list_samples.size() <= 7) {
        VENN_DIAGRAM(files, GENE_NAMES.out.collect())
    }

    if (params.concat) {
        FINAL_REPORT(CONCATENATE_BAMS.out.collect(), METHYLDACKEL.out.collect(), CORRELATION_AND_CLUSTERING.out.collect(), METILENE.out.collect(), params.cutoff_heatmap, samples_names)
    }
    else {
        FINAL_REPORT([], METHYLDACKEL.out.collect(), CORRELATION_AND_CLUSTERING.out.collect(), METILENE.out.collect(), params.cutoff_heatmap, samples_names)
    }
}
