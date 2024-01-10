# DNA-methylation-analysis

This RRBS Nextflow pipeline was created to discover the genes associated with differentially methylated regions from the CpG methylation patterns using [MethylDackel](https://github.com/dpryan79/MethylDackel) and [Metilene](http://legacy.bioinf.uni-leipzig.de/Software/metilene/).

The pipeline inputs **BAM files**, and outputs multiple txt and bedGraph files (according to the number of samples):
- Per base methylation metrics (.bedGraph)
- Differentially methylated regions (.bedGraph)
- Correlation matrix and PCA (.png)
- Heatmap with signature differences between the controls and samples (.pdf)
- Genomic distribution of CpGs with different methylation frequencies between samples and controls (.png)
- Genomic distribution of differentially methylated regions (.png)
- Closest RefSeq genes **(version from 2023-11-24)** to the differentially methylated regions (.txt/.bedGraph)
- Venn diagram of the closest genes (only if two or more samples were used as input) (.png)

![image](https://i.ibb.co/80RgK03/test.png)

# Install conda environment

To use this pipeline you need to have installed [conda](https://docs.conda.io/projects/conda/en/stable/user-guide/install/linux.html) and [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html).

```shell
git clone https://github.com/AnaValente/DNA-methylation-analysis/
cd DNA-methylation-analysis
conda env create -f methylation_env.yml
conda activate methylation
```

# Usage

Mandatory inputs:
 - `--files` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Path to scripts and samples folder
 - `--samples` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [String] sample names separated by comma (**always write the control name first!**)
 - `--replicates` &nbsp; [Integer] number of sample replicates

Optional inputs:
- `--concat ` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Option to concatenate BAM files from different runs
- `--cell_tpm` &nbsp;&nbsp;&nbsp; Optional file containing two collumns, one with gene names and other with expression levels in transcripts per million (TPM) specific to cell type (available in: https://www.ebi.ac.uk/gxa/experiments/E-MTAB-2770/Results) for gene name filtering
- `--cutoff_regions` &nbsp;&nbsp; [Integer] cutoff (from 1 to 100) for the difference between samples methylation frequency vs control methylation frequency for genomic annotations (default: 75)
- `--cutoff_heatmap` &nbsp;&nbsp; [Integer] cutoff (from 1 to 100) for the difference between samples methylation frenquecy vs control methylation frequency for clustering analysis (default: 100)

### Examples

**Example**
```
nextflow run Methylation_pipeline.nf --files /Scripts --samples 'Control','Sample1','Sample2'
```

**Example with BAM concatenation and genes filtered by file**
```
nextflow run Methylation_pipeline.nf --files /Scripts --samples 'Control','Sample1','Sample2' --cell_tpm E-MTAB-2770-query-results.tsv --concat 
```

**Example with different cutoffs**
```
nextflow run Methylation_pipeline.nf --files /Scripts --samples 'Control','Sample1','Sample2' --cutoff_regions 50 --cutoff_heatmap 75
```









