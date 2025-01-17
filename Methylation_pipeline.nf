#!/usr/bin/env nextflow

//params
params.files = false
params.samples = false
params.concat = false
params.replicates = false
params.cellTpm = false
params.cutoffRegions = 75
params.cutoffHeatmap = 100
params.genome = false
params.refseqGenes = false
params.noHg38 = false
params.inclusionBed = false
params.coverageFilter = 0
params.help = false

c_green = "\033[0;38;5;28m"
c_yellow = "\033[0;38;5;142m"
c_white = "\033[38;5;251m"
c_orange = "\033[38;5;202m"
c_blue = "\033[0;38;5;33m"
c_reset = "\033[0m"

log.info"""
        . ,-"-.   ,-"-. ,-"-.   ,-"-. ,-"-.   ,-"-. ,-"-.    ,-"-. ,-"-.   ,
         / ${c_green}| ${c_yellow}| ${c_reset}\\ / ${c_orange}| ${c_blue}| ${c_reset}/ ${c_blue}| ${c_green}| ${c_reset}\\ / ${c_yellow}| ${c_blue}| ${c_reset}/ ${c_yellow}| ${c_orange}| ${c_reset}\\ / ${c_green}| ${c_blue}| ${c_reset}/ ${c_orange}| ${c_blue}| ${c_reset}\\  / ${c_yellow}| ${c_green}| ${c_reset}/ ${c_orange}| ${c_blue}| ${c_reset}\\ /${c_reset}
          D N A  M E T H Y L A T I O N  A N A L Y S I S  P I P E L I N E
        / \\${c_blue}| ${c_orange}| ${c_reset}|\\| ${c_yellow}| ${c_green}|${c_reset}/ \\${c_green}| ${c_blue}| ${c_reset}|\\| ${c_orange}| ${c_green}|${c_reset}/ \\${c_orange}| ${c_yellow}| ${c_reset}|\\| ${c_blue}| ${c_green}|${c_reset}/ \\${c_yellow}| ${c_green}| ${c_reset}|\\|  ${c_orange}| ${c_blue}|${c_reset}/ \\${c_yellow}| ${c_green}| ${c_reset}|\\${c_reset}
           `-!-' `-!-'   `-!-' `-!-'   `-!-' `-!-'   `-!-'  `-!-'   `-!-' `-
        """.stripIndent()

if( params.help ) {

help = """
Usage:
    nextflow run Methylation_pipeline.nf --files <path> --samples <string> --replicates <integer> --genome <path> --refseqGenes <path> [options]

Mandatory arguments:
 --files 
    Path to scripts and samples folder.
 --samples 
    [String] sample names separated by comma (always write the control name first!).
 --replicates
    [Integer] number of sample replicates.
 --genome
    Path to the hg38 or other reference genome file (.fa.gz) (available in: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz).
 --refseqGenes
    Path to the refseq genes file in BED format (version from 2024-09-11 can be found in scripts folder for hg38 genome in GitHub).

Optional arguments:
 --concat
    Option to concatenate BAM files from different runs.
 --cellTpm
    Optional file containing two collumns, one with gene names and the other with expression levels in transcripts per million (TPM) for a cell line or cell type identical or similar to the cells under study (available in: https://www.ebi.ac.uk/gxa/experiments/E-MTAB-2770/Results) for gene name filtering.
 --noHg38
    Option to skip genomic regions output if not using hg38 reference genome.
 --cutoffRegions
    [Integer] cutoff (from 1 to 100) for the difference between samples methylation frequency vs control methylation frequency for genomic annotations [default: ${params.cutoffRegions}].
 --cutoffHeatmap
    [Integer] cutoff (from 1 to 100) for the difference between samples methylation frequency vs control methylation frequency for clustering analysis [default: ${params.cutoffHeatmap}].
 --inclusionBed
    Optional BED file containing regions of interest (A BED file containing CpG islands regions for hg38 reference genome can be found in scripts folder).

!!IMPORTANT!!
> All samples and additional files must be placed in the scripts folder.
> BAM files should start with the samples name provided in the –-samples input, followed by replicate number.
> If there are no replicates, input the number 1.
> Example:
    File names: Control_rep1_library.bam; Sample1_rep1_library.bam
    Pipeline: --samples "Control","Sample1" --replicates 1
""".stripMargin()

println(help)
exit(0)
}

process REFERENCE_GENOME {
    // Download fasta file of the reference genome 
    input:
    path genome
    
    output:
    path "*"

    """
    gunzip -c ${genome} > genome.fa
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
    path "${samples}*__library_sorted.b*"

    shell:
    '''
    for n in $(seq 1 !{n_replicates})
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
    path genome
    val samples
    val n_replicates
    val inclusion_bed

    output:
    path "*.bedGraph"

    publishDir './Pipeline_results/Methyldackel', pattern:"*svg" , mode: 'copy'

    """
    ./methyldackel.sh ${n_replicates} ${samples} ${inclusion_bed}
    """
}

process CORRELATION_AND_CLUSTERING {
    // correlation and clustering analysis

	publishDir './Pipeline_results', pattern:"*png" , mode: 'copy'
    publishDir './Pipeline_results/Methyldackel', pattern:"*CpGs.bedGraph" , mode: 'copy'

    input:
    path folder
    path bedgraphs
    val samples
    val number_samples
    val control
    val cutoff_regions
    val cutoff_heatmap
    val n_replicates
    val coverage_filter

    output:
    path "*"

    """
	python3 methyldackel_cluster_table.py ${bedgraphs} ${samples} ${number_samples} ${control} ${cutoff_regions} ${cutoff_heatmap} ${n_replicates} ${coverage_filter} && ./correlation_analysis.R correlation.csv ${number_samples} ${n_replicates} ${control} && ./heatmap.R diff_methylation_filtered.csv 
    """
}

process METILENE { 
    // Calculation of the differentially methylated regions
	
    input:
    path folder
    path bedgraph
    val control
    val samples
    val n_replicates

    output:
    path "*_qval.0.05.bedgraph"

    shell:
    '''
    samples_rep=""
    control_rep=""

    if  [[ ! -f !{control}_rep2_library_trimmed_CpGs.bedGraph ]] && [[ -f !{samples}_rep2_library_trimmed_CpGs.bedGraph ]]; then
        metilene_input.pl --in1 !{samples}_rep1_library_trimmed_CpGs.bedGraph,!{samples}_rep2_library_trimmed_CpGs.bedGraph --in2 !{control}_rep1_library_trimmed_CpGs.bedGraph
        metilene -a g1 -b g2 metilene_g1_g2.input | sort -V -k1,1 -k2,2n > metilene_!{samples}.bedgraph 
        metilene_output.pl -q metilene_!{samples}.bedgraph -o ./metilene_!{samples}_

    else

        for i in $(seq 1 !{n_replicates})
        do
            samples_rep=${samples_rep}",!{samples}_rep${i}_library_trimmed_CpGs.bedGraph"
            control_rep=${control_rep}",!{control}_rep${i}_library_trimmed_CpGs.bedGraph"
        done

        metilene_input.pl --in1 ${samples_rep/,/} --in2 ${control_rep/,/}
        metilene -a g1 -b g2 metilene_g1_g2.input | sort -V -k1,1 -k2,2n > metilene_!{samples}.bedgraph 
        metilene_output.pl -q metilene_!{samples}.bedgraph -o ./metilene_!{samples}_

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
    
    input:
    path folder
    path bedgraph
    val control
    val samples
    val cell_tpm
    val refseq_genes
    

    output:
    path "*_merged_genes*.*"

    """
	./get_refseq_gene_names.sh metilene_${samples}__qval.0.05.bedgraph ${refseq_genes} ${samples} && python3 expressed_cell_genes.py ${samples} ${cell_tpm}
    """
}

process VENN_DIAGRAM {
    // Venn Diagram

    publishDir './Pipeline_results', mode: 'copy'
    
    input:
    path folder
    path genes


    output:
    path "*.png"

    """
	./venn_diagram.R
    """
}

process FINAL_REPORT {
    
    publishDir './Pipeline_results/', pattern:"*_report.txt", mode: 'copy'
    publishDir './Pipeline_results/Metilene', pattern:"*_qval.0.05.bedgraph" , mode: 'copy'
    publishDir './Pipeline_results/Gene_names', pattern:"*_merged_genes*.*" , mode: 'copy'

    input:
    path bams
    path methyldackel
    path metilene
    path correlation
    path genes
    val cutoff_heatmap
    val samples

    output:
    path "*"

    shell:
    '''
    if compgen -G "*.bam" > /dev/null; then
        echo -e "Concatenated BAMS:\\n==================" >> Final_report.txt
        for file in *.bam; do
            echo -e "${file}:" >> Final_report.txt
            samtools flagstat ${file} >> Final_report.txt
        done
    fi

    echo -e "Methyldackel:\\n=============" >> Final_report.txt
    wc -l *_library_sorted_CpG.bedGraph | sed 's/_library.*//g' | awk '{print $2 ":", $1-1, "CpG sites\\n"}' >> Final_report.txt
    echo -e "\\nCorrelation and clustering:\\n===========================" >> Final_report.txt
    echo -e "$(wc -l correlation.csv | awk '{print $1-1}')" "CpG sites were used for the correlation matrix and" "$(python -c "print(round($(wc -l correlation.csv | awk '{print $1-1}') * 0.1))")" "CpG sites the for Principal Component Analysis\\n" >> Final_report.txt
    echo -e "$(wc -l diff_methylation_filtered.csv | awk '{print $1-1}')" "CpG sites were used for the hierarchical clustering with a difference of !{cutoff_heatmap}% between the controls and samples\\n" >> Final_report.txt
    echo -e "\\nMetilene:\\n=========" >> Final_report.txt
    
    for sample in !{samples}; do
        echo -ne "$(wc -l metilene_${sample}*_qval.0.05.bedgraph | sed 's/metilene_//' | sed 's/__.*//g' | awk '{print $2 ":", $1, "differentially methylated regions,"}') where $(awk '$4 > 0 {print ;}' metilene_${sample}*_qval.0.05.bedgraph  | wc -l) are hypermetilated" >> Final_report.txt
        echo -ne " and" "$(awk '$4 < 0 {print ;}' metilene_${sample}*_qval.0.05.bedgraph  | wc -l)" "are hypometilated\\n" >> Final_report.txt
    done

    echo -e "\\nGenes:\\n======" >> Final_report.txt

    for sample in !{samples}; do
        wc -l ${sample}*.txt | sed 's/__.*//g' | awk '{print $2 ":", $1, "genes"}' >> Final_report.txt
    done

    for i in *__* ; do cp "${i}" "${i/__/_}" ; done
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
        error("No reference genome file provided --genome")
    }

    if (!params.replicates) {
        error("No replicate number provided, if there are no replicates input 1 --replicates")
    }
    
    if (!params.refseqGenes && params.refseqGenes < 0) {
        error("No RefSeq file provided --refseqGenes")
    }
    
    if (params.cutoffRegions < 0 || params.cutoffRegions > 100 || params.cutoffHeatmap < 0 || params.cutoffHeatmap > 100) {
        error("Please provide a cutoff number between 0 and 100")
    }

    files = Channel.fromPath(params.files).collect()
    list_samples = params.samples.split(',') as List
    samples_chn = Channel.from(list_samples)
    control = list_samples.remove(0)
    samples_chn2 = Channel.from(list_samples)
    samples_names = list_samples.toString().replace("[", "").replace("]", "").replace(",", "")

    REFERENCE_GENOME(Channel.fromPath(params.genome))
    
    if (params.concat) {
        CONCATENATE_BAMS(files, samples_chn, params.replicates)
        METHYLDACKEL(files, CONCATENATE_BAMS.out.collect(), REFERENCE_GENOME.out.collect(), samples_chn, params.replicates, params.inclusionBed)
    } else {
        METHYLDACKEL(files, [], REFERENCE_GENOME.out.collect(), samples_chn, params.replicates, params.inclusionBed)
    }

    CORRELATION_AND_CLUSTERING(files, METHYLDACKEL.out.collect(), samples_names, list_samples.size(), control, params.cutoffRegions, params.cutoffHeatmap, params.replicates,params.coverageFilter)
    METILENE(files, CORRELATION_AND_CLUSTERING.out.collect(), control, samples_chn2, params.replicates)
    
    if (params.noHg38) {
        GENE_NAMES(files, METILENE.out.collect(), control, samples_chn2, params.cellTpm, params.refseqGenes)
    } else {
        GENOMIC_REGIONS(files, CORRELATION_AND_CLUSTERING.out.collect(), METILENE.out.collect())
        GENE_NAMES(files, METILENE.out.collect(), control, samples_chn2, params.cellTpm, params.refseqGenes)
    }

    if (list_samples.size() > 1 & list_samples.size() <= 6) {
        VENN_DIAGRAM(files, GENE_NAMES.out.collect())
    }

    if (params.concat) {
        FINAL_REPORT(CONCATENATE_BAMS.out.collect(), METHYLDACKEL.out.collect(), CORRELATION_AND_CLUSTERING.out.collect(), METILENE.out.collect(), GENE_NAMES.out.collect(), params.cutoffHeatmap, samples_names)
    }
    else {
        FINAL_REPORT([], METHYLDACKEL.out.collect(), CORRELATION_AND_CLUSTERING.out.collect(), METILENE.out.collect(), GENE_NAMES.out.collect(), params.cutoffHeatmap, samples_names)
    }
}