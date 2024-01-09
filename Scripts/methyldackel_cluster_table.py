#!/usr/bin/env python3

import sys
import glob
import pandas as pd

def trim_bedgraph(cutoff):

    dataframes = []
    for file in glob.glob("*CpG.bedGraph"):
        bedgraph = pd.read_table(file,skiprows=1, header=None, sep='\t')

        bedgraph = bedgraph.rename(columns={0:'chr',1:'start',2:'end',3:file.split('_library', 1)[0],4:'methylated_reads',5:'unmethylated_reads'})

        bedgraph['sum_reads'] = bedgraph['methylated_reads'] + bedgraph['unmethylated_reads']

        bedgraph_sum = bedgraph[bedgraph['sum_reads'] >= cutoff]

        bedgraph_sum = bedgraph_sum.drop(['methylated_reads','unmethylated_reads','sum_reads'], axis=1)

        dataframes.append(bedgraph_sum)

        merged_dataframe = dataframes[0]

        for df in dataframes[1:]:
            merged_dataframe = pd.merge(merged_dataframe, df, on=['chr','start','end'])

            merged_dataframe = pd.concat([merged_dataframe[['chr','start','end']],merged_dataframe[merged_dataframe.columns.difference(['chr','start','end'])].sort_index(axis=1)],ignore_index=False,axis=1)

    return merged_dataframe

def merge_n_cutoff(df,samples,n_samples,control,cutoff):

    included = pd.DataFrame()
    filtered = pd.DataFrame()

    samples_names = samples[-(int(n_samples)+4):-4]

    for i in range(len(df)):
        count = 0

        for j in samples_names:
            if control+'rep2' in df.columns:
                if (abs(int(df[j+'_rep1'].iloc[i]) - int(df[control+'_rep1'].iloc[i])) >= cutoff and abs(int(df[j+'_rep2'].iloc[i]) - int(df[control+'_rep1'].iloc[i])) >= cutoff and abs(int(df[j+'_rep1'].iloc[i]) - int(df[control+'_rep2'].iloc[i])) >= cutoff and abs(int(df[j+'_rep2'].iloc[i]) - int(df[control+'_rep2'].iloc[i])) >= cutoff):
                    included = pd.concat([included, df.iloc[[i]]])
                    count = count + 1
            elif (abs(int(df[j+'_rep1'].iloc[i]) - int(df[control+'_rep1'].iloc[i])) >= cutoff and abs(int(df[j+'_rep2'].iloc[i]) - int(df[control+'_rep1'].iloc[i])) >= cutoff):
                included = pd.concat([included, df.iloc[[i]]])
                count = count + 1

        if count == 1:
            filtered = pd.concat([filtered, df.iloc[[i]]])

    included = included.drop_duplicates().reset_index(drop=True)
    filtered = filtered.drop_duplicates().reset_index(drop=True)

    return included,filtered

def cutoff_genm_regions(df,samples,n_samples,control,cutoff):

    samples_names = samples[-(int(n_samples)+4):-4]

    for j in samples_names:
        nanomaterial = pd.DataFrame()

        for i in range(len(df)):
            if control+'rep2' in df.columns:
                if (abs(int(df[j+'_rep1'].iloc[i]) - int(df[control+'_rep1'].iloc[i])) >= cutoff and (abs(int(df[j+'_rep2'].iloc[i]) - int(df[control+'_rep1'].iloc[i]))) >= cutoff and abs(int(df[j+'_rep1'].iloc[i]) - int(df[control+'_rep2'].iloc[i])) >= cutoff and (abs(int(df[j+'_rep2'].iloc[i]) - int(df[control+'_rep2'].iloc[i]))) >= cutoff):
                    nanomaterial = pd.concat([nanomaterial, df.iloc[[i]]])
                    nanomaterial.to_csv(j+'_methylation.csv', sep='\t', index=False, header=False)

            elif (abs(int(df[j+'_rep1'].iloc[i]) - int(df[control+'_rep1'].iloc[i])) >= cutoff and (abs(int(df[j+'_rep2'].iloc[i]) - int(df[control+'_rep1'].iloc[i]))) >= cutoff):
                nanomaterial = pd.concat([nanomaterial, df.iloc[[i]]])
                nanomaterial.to_csv(j+'_methylation.csv', sep='\t', index=False, header=False)


if __name__ == "__main__":

    correlation_df = trim_bedgraph(0)
    correlation_df.to_csv('correlation.csv', sep='\t', index=False)

    included_df,unfiltered_df = merge_n_cutoff(trim_bedgraph(5),list(sys.argv),sys.argv[len(sys.argv)-4],sys.argv[len(sys.argv)-3],int(sys.argv[len(sys.argv)-2]))
    nocluded_df,filtered_df = merge_n_cutoff(trim_bedgraph(5),list(sys.argv),sys.argv[len(sys.argv)-4],sys.argv[len(sys.argv)-3],int(sys.argv[len(sys.argv)-1]))

    filtered_df.to_csv('diff_methylation_filtered.csv', sep='\t', index=False, header=True)
    included_df.to_csv('diff_methylation_all.csv', sep='\t', index=False, header=True)

    cutoff_genm_regions(included_df,list(sys.argv), sys.argv[len(sys.argv)-4],sys.argv[len(sys.argv)-3],int(sys.argv[len(sys.argv)-2]))
