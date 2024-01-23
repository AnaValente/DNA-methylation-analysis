#!/usr/bin/env python3

# imports
import sys
from gprofiler import GProfiler
import pandas as pd

def cell_common_genes(df,tpm,tpm5):

    common_genes = pd.DataFrame()
    excluded_genes = pd.DataFrame()

    for i in range(len(df)):
        if df['mean_methylation_dif'][i] < 0:

            if tpm.isin([df['Gene Name'][i]]).any().any() == True:
                common_genes = pd.concat([common_genes, df.iloc[[i]]]) # store the common genes between files in dataframe
            else:
                excluded_genes = pd.concat([excluded_genes, df.iloc[[i]]])
        else:
            if tpm5.isin([df['Gene Name'][i]]).any().any() == True:

                common_genes = pd.concat([common_genes, df.iloc[[i]]]) # store the common genes between files in dataframe
            else:
                excluded_genes = pd.concat([excluded_genes, df.iloc[[i]]])

    return common_genes.reset_index(drop=True)

if __name__ == "__main__":

    if sys.argv[3] != 'false':
        cell_tpm = pd.read_csv(sys.argv[3], skiprows=4, sep='\t')
        
        # convert the df gene names to HGNC nomenclature
        gp = GProfiler(return_dataframe=True)
        gp = gp.convert(organism='hsapiens',query=cell_tpm['Gene ID'].values.tolist(),target_namespace='HGNC')
        cell_tpm = cell_tpm.join(gp['converted']).drop(['Gene Name'],axis=1)
        cell_tpm5 = cell_tpm[cell_tpm.iloc[:, 1] >= 5.0] # select only the genes with more than 5 TPM

        ctrl_NM_genes = pd.read_csv(sys.argv[2]+'_'+sys.argv[1]+'_refseq.bedgraph', sep='\t',header=None,names=['chr','start','end','mean_methylation_dif','Gene Name','dist'])
        common_genes_NM = cell_common_genes(ctrl_NM_genes,cell_tpm,cell_tpm5)

        common_genes_NM.to_csv(sys.argv[2]+'_'+sys.argv[1]+'_merged_genes.bedgraph',index=False,header=False,sep='\t')
        common_genes_NM['Gene Name'].to_csv(sys.argv[2]+'_'+sys.argv[1]+'_merged_genes.txt',index=False,header=False,sep='\t')
    else:
        ctrl_NM_genes = pd.read_csv(sys.argv[2]+'_'+sys.argv[1]+'_refseq.bedgraph', sep='\t',header=None,names=['chr','start','end','mean_methylation_dif','Gene Name','dist'])
        ctrl_NM_genes.to_csv(sys.argv[2]+'_'+sys.argv[1]+'_merged_genes.bedgraph',index=False,header=False,sep='\t')
        ctrl_NM_genes['Gene Name'].to_csv(sys.argv[2]+'_'+sys.argv[1]+'_merged_genes.txt',index=False,header=False,sep='\t')
