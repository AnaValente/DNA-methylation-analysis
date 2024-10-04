#!/usr/bin/env python3

# imports
import sys
from gprofiler import GProfiler
import pandas as pd

def cell_common_genes(df, tpm, tpm5):
    """
    Find common genes based on methylation frequency and TPM values.
    """
    common_genes = pd.DataFrame()  # DataFrame to store common genes
    excluded_genes = pd.DataFrame()  # DataFrame to store excluded genes

    for i in range(len(df)):
        if df['mean_methylation_dif'][i] < 0:

            if tpm.isin([df['Gene Name'][i]]).any().any() == True:
                common_genes = pd.concat([common_genes, df.iloc[[i]]])  # Store common genes in dataframe
            else:
                excluded_genes = pd.concat([excluded_genes, df.iloc[[i]]])  # Store excluded genes
        else:
            if tpm5.isin([df['Gene Name'][i]]).any().any() == True:
                common_genes = pd.concat([common_genes, df.iloc[[i]]])  
            else:
                excluded_genes = pd.concat([excluded_genes, df.iloc[[i]]]) 

    return common_genes.reset_index(drop=True), excluded_genes.reset_index(drop=True)


if __name__ == "__main__":

    if sys.argv[2] != 'false':
        cell_tpm = pd.read_csv(sys.argv[2], skiprows=4, sep='\t') # Read TPM table
        
        # Convert the gene names to HGNC nomenclature
        gp = GProfiler(return_dataframe=True)
        gp = gp.convert(organism='hsapiens', query=cell_tpm['Gene ID'].values.tolist(), target_namespace='HGNC')
        cell_tpm = cell_tpm.join(gp['converted']).drop(['Gene Name'], axis=1)
        cell_tpm5 = cell_tpm[cell_tpm.iloc[:, 1] >= 5.0]  # Select genes with TPM >= 5

        ctrl_NM_genes = pd.read_csv(sys.argv[1]+'_refseq.bedgraph', sep='\t', header=None, names=['chr','start','end','mean_methylation_dif','Gene Name','dist'])
        common_genes_NM,excluded_genes = cell_common_genes(ctrl_NM_genes, cell_tpm, cell_tpm5)

        excluded_genes.to_csv(sys.argv[1]+'__merged_genes_excluded.bedgraph', index=False, header=False, sep='\t')
        common_genes_NM.to_csv(sys.argv[1]+'__merged_genes.bedgraph', index=False, header=False, sep='\t')
        common_genes_NM['Gene Name'].to_csv(sys.argv[1]+'__merged_genes.txt', index=False, header=False, sep='\t')
    else:

        ctrl_NM_genes = pd.read_csv(sys.argv[1]+'__refseq.bedgraph', sep='\t', header=None, names=['chr','start','end','mean_methylation_dif','Gene Name','dist'])
        ctrl_NM_genes.to_csv(sys.argv[1]+'__merged_genes.bedgraph', index=False, header=False, sep='\t')
        ctrl_NM_genes['Gene Name'].to_csv(sys.argv[1]+'__merged_genes.txt', index=False, header=False, sep='\t')
