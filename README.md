# DNA-methylation-analysis

This RRBS Nextflow pipeline was created to discover the differentialy methylated genes from the CpG methylation patterns using [MethylDackel](https://github.com/dpryan79/MethylDackel) and [Metilene](http://legacy.bioinf.uni-leipzig.de/Software/metilene/).

The pipeline inputs **BAM files of two replicates** (from multiple runs if needed), and outputs multiple txt and bedGraph files (according to the number of samples) with the gene names of the different methylation regions (DMR), as well as:
- MethylDackel per base methylation metrics (.bedGraph)
- Metilene different methylated regions (.bedGraph)
- Replicates correlation matrix and PCA (.png)
- Heatmap with "signature" differences between the control and samples (.pdf)
- Genomic regions distribution of the methylation metrics and DMRs (.png)
- Genes venn diagram (.png)

![image](https://i.ibb.co/p1zYpTc/test.png)

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
 - `--files` &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Path to scripts and samples folder
 - `--samples` &nbsp;&nbsp; [String] sample names separated by comma (**always write the control name first!**)
 - `--cell_tpm` &nbsp; Path to file containing cell gene expresion levels in transcripts per million (TPM) (tsv format)

Optional inputs:
- `--concat` &nbsp;&nbsp;&nbsp; Option to concatenate BAM files from different runs
- `--cutoff_regions` &nbsp;&nbsp; [Integer] cutoff (from 0 to 100) for the difference between samples methylation metrics vs control methylation metrics for genomic annotations (default: 75)
- `--cutoff_heatmap` &nbsp;&nbsp; [Integer] cutoff (from 0 to 100) for the difference between samples methylation metrics vs control methylation metrics for clustering analysis (default: 100)

### Examples

**Example**
```
nextflow run Methylation_pipeline.nf --files /Scripts --samples 'Control','Sample1','Sample2','Sample3' --cell_tpm E-MTAB-2770-query-results.tsv
```

**Example with BAM concatenation**
```
nextflow run Methylation_pipeline.nf --files /Scripts --samples 'Control','Sample1','Sample2','Sample3' --cell_tpm E-MTAB-2770-query-results.tsv --concat
```

**Example with different cutoffs**
```
nextflow run Methylation_pipeline.nf --files /Scripts --samples 'Control','Sample1','Sample2','Sample3' --cell_tpm E-MTAB-2770-query-results.tsv --cutoff_regions 50 --cutoff_heatmap 75
```









