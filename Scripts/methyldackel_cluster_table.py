#!/usr/bin/env python3

import sys
import glob
import pandas as pd
from statistics import mean

def trim_bedgraph(cutoff,control):

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

    control_cols = [col for col in merged_dataframe.columns if control in col]
    first_cols = ['chr','start','end'] + control_cols

    merged_dataframe = pd.concat([merged_dataframe[first_cols],merged_dataframe[merged_dataframe.columns.difference(first_cols)].sort_index(axis=1)],ignore_index=False,axis=1)
    
    return merged_dataframe

def merge_n_cutoff(df,n_replicates,samples,n_samples,control,cutoff):
    included = pd.DataFrame()
    filtered = pd.DataFrame()

    samples_names = samples[-(int(n_samples)+5):-5]

    for i in range(len(df)):
        count = 0

        for j in samples_names:

            diff_samples = list()
            diff_control = list()

            if j+'_rep2' in df.columns and not control+'_rep2' in df.columns:
                if (abs(int(df[j+'_rep1'].iloc[i]) - int(df[control+'_rep1'].iloc[i])) >= cutoff and abs(int(df[j+'_rep2'].iloc[i]) - int(df[control+'_rep1'].iloc[i])) >= cutoff):
                    included = pd.concat([included, df.iloc[[i]]])
                    count = count + 1
            else:
                for r in range(1,n_replicates+1):
                    diff_samples.append(int(df[j + "_rep" + str(r)].iloc[i]))
                    diff_control.append(int(df[control + "_rep" + str(r)].iloc[i]))
                
                if abs(mean(diff_samples) - mean(diff_control)) >= cutoff:
                    included = pd.concat([included, df.iloc[[i]]])
                    count = count + 1

        if count == 1:
            filtered = pd.concat([filtered, df.iloc[[i]]])

    included = included.drop_duplicates().reset_index(drop=True)
    filtered = filtered.drop_duplicates().reset_index(drop=True)

    return included,filtered

def cutoff_genm_regions(df,n_replicates,samples,n_samples,control,cutoff):

    samples_names = samples[-(int(n_samples)+5):-5]

    for j in samples_names:
        
        nanomaterial = pd.DataFrame()

        for i in range(len(df)):

            diff_samples = list()
            diff_control = list()

            if j+'_rep2' in df.columns and not control+'_rep2' in df.columns:
                if (abs(int(df[j+'_rep1'].iloc[i]) - int(df[control+'_rep1'].iloc[i])) >= cutoff and abs(int(df[j+'_rep2'].iloc[i]) - int(df[control+'_rep1'].iloc[i])) >= cutoff):
                    nanomaterial = pd.concat([nanomaterial, df.iloc[[i]]])
            
            else:
                for r in range(1,n_replicates+1):
                    diff_samples.append(int(df[j+"_rep" + str(r)].iloc[i]))
                    diff_control.append(int(df[control+"_rep" + str(r)].iloc[i]))

                if abs(mean(diff_samples) - mean(diff_control)) >= cutoff:
                    nanomaterial = pd.concat([nanomaterial, df.iloc[[i]]])
        
        nanomaterial.to_csv(j+'_methylation.csv', sep='\t', index=False, header=False)

if __name__ == "__main__":

    correlation_df = trim_bedgraph(0, sys.argv[len(sys.argv)-4])
    correlation_df.to_csv('correlation.csv', sep='\t', index=False)

    trim_bedgraph_5 = trim_bedgraph(5, sys.argv[len(sys.argv)-4])

    included_df,unfiltered_df = merge_n_cutoff(trim_bedgraph_5,int(sys.argv[len(sys.argv)-1]),list(sys.argv),sys.argv[len(sys.argv)-5],sys.argv[len(sys.argv)-4],int(sys.argv[len(sys.argv)-3]))
    removed_df,filtered_df = merge_n_cutoff(trim_bedgraph_5,int(sys.argv[len(sys.argv)-1]),list(sys.argv),sys.argv[len(sys.argv)-5],sys.argv[len(sys.argv)-4],int(sys.argv[len(sys.argv)-2]))

    filtered_df.to_csv('diff_methylation_filtered.csv', sep='\t', index=False, header=True)
    included_df.to_csv('diff_methylation_all.csv', sep='\t', index=False, header=True)

    cutoff_genm_regions(included_df,int(sys.argv[len(sys.argv)-1]),list(sys.argv),sys.argv[len(sys.argv)-5],sys.argv[len(sys.argv)-4],int(sys.argv[len(sys.argv)-3]))