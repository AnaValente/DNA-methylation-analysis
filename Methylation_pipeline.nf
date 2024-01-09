#!/usr/bin/env nextflow

log.info"""\

         M E T H Y L A T I O N   A N A L Y S I S   P I P E L I N E
         =========================================================
         """

//params
params.files = false
params.samples = false
params.concat = false
params.cell_tpm = 'E-MTAB-2770-query-results.tsv'
params.cutoff_regions = 75
params.cutoff_heatmap = 100

process REFERENCE_GENOME_HG38 {
    // Download fasta file of the reference genome 
    
    publishDir "./Pipeline_results/Methyldackel", mode: "copy"

    output:
    path "*"

    """
    curl http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | gunzip -c > hg38.fa
    """
}

process CONCATENATE_BAMS_REP1 {
    // Concatenate the replicate 1 bam files for each sample

    publishDir "./Pipeline_results/Concatenated_files", mode: "copy"
    

    input:
    path folder
    val samples

    output:
    path "${samples}*_library_sorted.b*"

    """
    samtools merge ${samples}_rep1_library.bam ${samples}_rep1*.bam && samtools sort ${samples}_rep1_library.bam -o ${samples}_rep1_library_sorted.bam -O BAM && samtools index ${samples}_rep1_library_sorted.bam
    samtools flagstat ${samples}_rep1_library_sorted.bam > ${samples}_rep1_report.txt
    """
}

process CONCATENATE_BAMS_REP2 {
    // Concatenate the replicate 2 bam files for each sample

    publishDir "./Pipeline_results/Concatenated_files", mode: "copy"
    

    input:
    path folder
    val nano

    output:
    path "${nano}*_library_sorted.b*"
    
    """
    samtools merge ${nano}_rep2_library.bam ${nano}_rep2*.bam && samtools sort ${nano}_rep2_library.bam -o ${nano}_rep2_library_sorted.bam -O BAM && samtools index ${nano}_rep2_library_sorted.bam
    samtools flagstat ${nano}_rep2_library_sorted.bam > ${nano}_rep2_report.txt
    """
}

process METHYLDACKEL_REP1 {
    // Methylation calling for replicate 1

    publishDir './Pipeline_results/Methyldackel', mode: 'copy'
    
    input:
    path folder
    path hg38
    path bam_rep1
    val samples


    output:
    path "*.bedGraph"

    """
    if [ -f ${samples}_rep1_library_sorted.bam ]; then
        MethylDackel extract --mergeContext -l CpGIslands_mod.bed hg38.fa ${samples}_rep1_library_sorted.bam
        MethylDackel extract --mergeContext -l CpGIslands_mod.bed --fraction hg38.fa ${samples}_rep1_library_sorted.bam
    else
        mv ${samples}_rep1*.bam ${samples}_rep1_library_sorted.bam
        MethylDackel extract --mergeContext -l CpGIslands_mod.bed hg38.fa ${samples}_rep1_library_sorted.bam
        MethylDackel extract --mergeContext -l CpGIslands_mod.bed --fraction hg38.fa ${samples}_rep1_library_sorted.bam
    fi    
    """
}

process METHYLDACKEL_REP2 {
    // Methylation calling for replicate 2

    publishDir './Pipeline_results/Methyldackel', mode: 'copy'
    
    input:
    path folder
    path hg38
    path bam_rep1
    val samples

    output:
    path "*.bedGraph"

    """
    if [ -f ${samples}_rep2_library_sorted.bam ]; then
        MethylDackel extract --mergeContext -l CpGIslands_mod.bed hg38.fa ${samples}_rep2_library_sorted.bam
        MethylDackel extract --mergeContext -l CpGIslands_mod.bed --fraction hg38.fa ${samples}_rep2_library_sorted.bam
    else
        mv ${samples}_rep2*.bam ${samples}_rep2_library_sorted.bam
        MethylDackel extract --mergeContext -l CpGIslands_mod.bed hg38.fa ${samples}_rep2_library_sorted.bam
        MethylDackel extract --mergeContext -l CpGIslands_mod.bed --fraction hg38.fa ${samples}_rep2_library_sorted.bam
    fi    
    """
}

process CORRELATION_AND_CLUSTERING {
    // correlation and clustering analysis

    publishDir './Pipeline_results', mode: 'copy'
	
    input:
    path folder
    path bedgraph_rep1
    path bedgraph_rep2
    val samples
    val number_samples
    val control
    val cutoff_regions
    val cutoff_heatmap

    output:
    path "*"

    """
	python3 methyldackel_cluster_table.py ${bedgraph_rep1} ${bedgraph_rep2} ${samples} ${number_samples} ${control} ${cutoff_regions} ${cutoff_heatmap} && ./correlation_analysis.R correlation.csv && ./heatmap.R diff_methylation_filtered.csv 
    """
}

process METILENE { 
    // Calculation of the differentially methylated regions

    publishDir './Pipeline_results/Metilene', mode: 'copy'
    
	
    input:
    path bedgraph_rep1
    path bedgraph_rep2
    val control
    val samples

    output:
    path "*_qval.0.05.bedgraph"

    """
    if [ -f ${control}_rep2*.bedGraph ]; then
        metilene_input.pl --in1 ${samples}_rep1_library_sorted_CpG.meth.bedGraph,${samples}_rep2_library_sorted_CpG.meth.bedGraph --in2 ${control}_rep1_library_sorted_CpG.meth.bedGraph,${control}_rep2_library_sorted_CpG.meth.bedGraph
        metilene -a g1 -b g2 metilene_g1_g2.input | sort -V -k1,1 -k2,2n > metilene_${samples}_${control}.bedgraph
        metilene_output.pl -q metilene_${samples}_${control}.bedgraph -o ./metilene_${samples}_${control}
    else
        metilene_input.pl --in1 ${samples}_rep1_library_sorted_CpG.meth.bedGraph,${samples}_rep2_library_sorted_CpG.meth.bedGraph --in2 ${control}_rep1_library_sorted_CpG.meth.bedGraph
        metilene -a g1 -b g2 metilene_g1_g2.input | sort -V -k1,1 -k2,2n > metilene_${samples}_${control}.bedgraph
        metilene_output.pl -q metilene_${samples}_${control}.bedgraph -o ./metilene_${samples}_${control}
    fi
    """
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

    publishDir './Pipeline_results/Gene_Names', mode: 'copy'
    
	
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

workflow {
    files = Channel.fromPath(params.files).collect()
    list_samples = params.samples.split(',') as List
    samples_chn = Channel.from(list_samples)
    control = list_samples.remove(0)
    samples_chn2 = Channel.from(list_samples)
    samples_names = list_samples.toString().replace("[", "").replace("]", "").replace(",", "")

    REFERENCE_GENOME_HG38()
    if (params.concat) {
        CONCATENATE_BAMS_REP1(files, samples_chn)
        CONCATENATE_BAMS_REP2(files, samples_chn)
        METHYLDACKEL_REP1(files,REFERENCE_GENOME_HG38.out.collect(),CONCATENATE_BAMS_REP1.out.collect(),samples_chn)
        METHYLDACKEL_REP2(files,REFERENCE_GENOME_HG38.out.collect(),CONCATENATE_BAMS_REP2.out.collect(),samples_chn)
        CORRELATION_AND_CLUSTERING(files, METHYLDACKEL_REP1.out.collect(), METHYLDACKEL_REP2.out.collect(),samples_names,list_samples.size(),control,params.cutoff_regions,params.cutoff_heatmap)
        METILENE(METHYLDACKEL_REP1.out.collect(), METHYLDACKEL_REP2.out.collect(), control, samples_chn2)
        GENOMIC_REGIONS(files, CORRELATION_AND_CLUSTERING.out.collect(), METILENE.out.collect())
        GENE_NAMES(files,METILENE.out.collect(), control, samples_chn2,params.cell_tpm)
        VENN_DIAGRAM(files, GENE_NAMES.out.collect())
        }
    else {
        METHYLDACKEL_REP1(files,REFERENCE_GENOME_HG38.out.collect(),[],samples_chn)
        METHYLDACKEL_REP2(files,REFERENCE_GENOME_HG38.out.collect(),[],samples_chn)
        CORRELATION_AND_CLUSTERING(files, METHYLDACKEL_REP1.out.collect(), METHYLDACKEL_REP2.out.collect(),samples_names,list_samples.size(),control,params.cutoff_regions,params.cutoff_heatmap)
        METILENE(METHYLDACKEL_REP1.out.collect(),METHYLDACKEL_REP2.out.collect(), control, samples_chn2)
        GENOMIC_REGIONS(files, CORRELATION_AND_CLUSTERING.out.collect(), METILENE.out.collect())
        GENE_NAMES(files, METILENE.out.collect(), control, samples_chn2,params.cell_tpm)
        VENN_DIAGRAM(files, GENE_NAMES.out.collect())
        }
}
