# DNA-methylation-analysis

This RRBS Nextflow pipeline was created to discover the genes associated with differentially methylated regions from the CpG methylation patterns using [MethylDackel](https://github.com/dpryan79/MethylDackel) and [Metilene](http://legacy.bioinf.uni-leipzig.de/Software/metilene/).

The pipeline inputs **BAM files**, and outputs multiple txt and bedGraph files (according to the number of samples):
- Per base methylation metrics (.bedGraph)
- Differentially methylated regions (.bedGraph)
- Correlation matrix and PCA (.png)
- Heatmap with signature differences between the controls and samples (.png)
- Genomic distribution across the **hg38 reference genome** of CpGs with different methylation frequencies between samples and controls (.png)
- Genomic distribution across the **hg38 reference genome** of differentially methylated regions (.png)
- Closest RefSeq genes **(version from 2023-11-24)** to the differentially methylated regions (.txt/.bedGraph)
- Venn diagram of the closest genes **(only if two or more samples were used as input)** (.png)

![image](https://i.ibb.co/80RgK03/test.png)

> [!NOTE]
> **This pipeline uses by default the hg38 reference genome**
>
> **Other reference genomes are available to used but the genomic regions output will be skipped**

# Install conda environment

To use this pipeline you need to have installed [conda](https://docs.conda.io/projects/conda/en/stable/user-guide/install/linux.html) and [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html).

```shell
git clone https://github.com/AnaValente/DNA-methylation-analysis/
cd DNA-methylation-analysis
conda env create -f methylation_env.yml
conda activate methylation
```

# Usage

**Mandatory inputs:**
 - `--files` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Path to scripts and samples folder
 - `--samples` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [String] sample names separated by comma **(always write the control name first!)**
 - `--replicates` &nbsp;&nbsp; [Integer] number of sample replicates
 - `--genome` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Path to the hg38 or other reference genome file **(.fa.gz)** (available in: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)
 - `--refseq_genes` Path to the refseq genes file in bed format **(version from 2023-11-24 can be found in scripts folder for hg38 genome)**

**Optional inputs:**
- `--concat ` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Option to concatenate BAM files from different runs
- `--cell_tpm` &nbsp;&nbsp;&nbsp; Optional file containing two collumns, one with gene names and the other with expression levels in transcripts per million (TPM) for a cell line or cell type identical or similar to the cells under study (available in: https://www.ebi.ac.uk/gxa/experiments/E-MTAB-2770/Results) for gene name filtering
-  `--no_hg38` &nbsp;&nbsp; Option to skip genomic regions output if not using hg38 reference genome
- `--cutoff_regions` &nbsp;&nbsp; [Integer] cutoff (from 1 to 100) for the difference between samples methylation frequency vs control methylation frequency for genomic annotations (default: 75)
- `--cutoff_heatmap` &nbsp;&nbsp; [Integer] cutoff (from 1 to 100) for the difference between samples methylation frequency vs control methylation frequency for clustering analysis (default: 100)

> [!IMPORTANT]
> **All samples and additional files must be placed in the scripts folder**
> 
> **BAM files should start with the samples name provided in the –-samples input, followed by replicate number**
>
>  **Example:**
>  - File names: **Control_rep1**_library.bam; **Sample1_rep1**_library.bam
>  - Pipeline : --samples **"Control","Sample1"** --replicates **1**

# Examples

**Example**
```
nextflow run Methylation_pipeline.nf --files "Scripts/*" --samples 'Control','Sample1','Sample2' --replicates 2 --genome hg38.fa.gz --refseq_genes RefSeq_genes.bed
```

**Example with BAM concatenation**
```
nextflow run Methylation_pipeline.nf --files "Scripts/*" --samples 'Control','Sample1','Sample2' --replicates 2 --genome hg38.fa.gz --refseq_genes RefSeq_genes.bed --concat
```

**Example with genes filtered by file**
```
nextflow run Methylation_pipeline.nf --files "Scripts/*" --samples 'Control','Sample1','Sample2' --replicates 2 --genome hg38.fa.gz --refseq_genes RefSeq_genes.bed --cell_tpm E-MTAB-2770-query-results.tsv 
```

**Example with different cutoffs**
```
nextflow run Methylation_pipeline.nf --files "Scripts/*" --samples 'Control','Sample1','Sample2' --replicates 2 --genome hg38.fa.gz --refseq_genes RefSeq_genes.bed --cutoff_regions 50 --cutoff_heatmap 75
```
