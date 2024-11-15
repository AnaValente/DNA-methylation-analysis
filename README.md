# DNA-methylation-analysis
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14169665.svg)](https://doi.org/10.5281/zenodo.14169665)

This RRBS Nextflow pipeline was created to discover the genes associated with differentially methylated regions from the CpG methylation patterns using [MethylDackel](https://github.com/dpryan79/MethylDackel) and [Metilene](http://legacy.bioinf.uni-leipzig.de/Software/metilene/).

The pipeline inputs **BAM files**, and outputs multiple txt and bedGraph files (according to the number of samples):
- Per base methylation metrics (.bedGraph)
- Differentially methylated regions (.bedGraph)
- Correlation matrix and PCA (.png)
- Heatmap with signature differences between the controls and samples (.png)
- Genomic distribution across the **hg38 reference genome** of CpGs with different methylation frequencies between samples and controls (.png)
- Genomic distribution across the **hg38 reference genome** of differentially methylated regions (.png)
- Closest RefSeq genes **(version from 2024-09-11)** to the differentially methylated regions (.txt/.bedGraph)
- Venn diagram of the closest genes **(only if two or more samples were used as input)** (.png)

![image](https://i.ibb.co/80RgK03/test.png)

> [!NOTE]
> **This pipeline uses by default the hg38 reference genome**
>
> **Other reference genomes are available to use, but the genomic regions output will be skipped**

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
 - `--refseqGenes` Path to the refseq genes file in BED format **(version from 2024-09-11 can be found in scripts folder for hg38 genome)**

**Optional inputs:**
- `--concat ` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Option to concatenate BAM files from different runs
- `--cellTpm` &nbsp;&nbsp;&nbsp; Optional file containing two collumns, one with gene names and the other with expression levels in transcripts per million (TPM) for a cell line or cell type identical or similar to the cells under study (available in: https://www.ebi.ac.uk/gxa/experiments/E-MTAB-2770/Results) for gene name filtering
-  `--noHg38` &nbsp;&nbsp; Option to skip genomic regions output if not using hg38 reference genome
- `--cutoffRegions` &nbsp;&nbsp; [Integer] cutoff (from 1 to 100) for the difference between samples methylation frequency vs control methylation frequency for genomic annotations (default: 75)
- `--cutoffHeatmap` &nbsp;&nbsp; [Integer] cutoff (from 1 to 100) for the difference between samples methylation frequency vs control methylation frequency for clustering analysis (default: 100)
- `--inclusionBed` &nbsp;&nbsp; Optional BED file containing regions of interest **(A BED file containing CpG islands regions for hg38 reference genome can be found in scripts folder)**

> [!IMPORTANT]
> - **All samples and additional files must be placed in the scripts folder**
> 
> - **BAM files should start with the samples name provided in the –-samples input, followed by replicate number**
>
> - **If there are no replicates, input the number 1**
>
>  - **Example:**
>    - **File names:** **Control_rep1**_library.bam; **Sample1_rep1**_library.bam
>    - **Pipeline :** --samples **"Control","Sample1"** --replicates **1**

# Examples

**Example**
```
nextflow run Methylation_pipeline.nf --files "Scripts/*" --samples 'Control','Sample1','Sample2' --replicates 2 --genome hg38.fa.gz --refseqGenes RefSeq_genes.bed
```

**Example with BAM concatenation**
```
nextflow run Methylation_pipeline.nf --files "Scripts/*" --samples 'Control','Sample1','Sample2' --replicates 2 --genome hg38.fa.gz --refseqGenes RefSeq_genes.bed --concat
```

**Example with genes filtered by file and regions of interest**
```
nextflow run Methylation_pipeline.nf --files "Scripts/*" --samples 'Control','Sample1','Sample2' --replicates 2 --genome hg38.fa.gz --refseqGenes RefSeq_genes.bed --cellTpm E-MTAB-2770-query-results.tsv --inclusionBed CpGIslands_mod.bed
```

**Example with different cutoffs**
```
nextflow run Methylation_pipeline.nf --files "Scripts/*" --samples 'Control','Sample1','Sample2' --replicates 2 --genome hg38.fa.gz --refseqGenes RefSeq_genes.bed --cutoffRegions 50 --cutoffHeatmap 75
```
