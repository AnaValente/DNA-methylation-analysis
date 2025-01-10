#!/usr/bin/env python3

import sys
import glob
import pandas as pd
from statistics import mean

def merge_dfs(df_list):
    ''''
    Merge list of Dataframes.
    '''
    merged_dfs = df_list[0]

    for df in df_list[1:]:
        merged_dfs = pd.merge(merged_dfs, df, on=['chr','start','end'], how='inner') # Merge DataFrames in List
    
    return merged_dfs

def trim_bedgraph(cutoff, control):
    ''''
    Merge and trim bedGraph files based on cutoff.
    '''
    dataframes = [] # List to store all DataFrames

    for file in glob.glob("*CpG.bedGraph"):
        prefix = file.split('_library', 1)[0]  
        bedgraph = pd.read_table(file, skiprows=1, header=None, sep='\t', names=['chr', 'start', 'end', prefix, 'methylated_reads', 'unmethylated_reads'])  # Set column names directly

        bedgraph['sum_reads'] = bedgraph['methylated_reads'] + bedgraph['unmethylated_reads']
        bedgraph_filtered = bedgraph[bedgraph['sum_reads'] >= cutoff] # Select number of reads based on the cutoff
        
        bedgraph_filtered = bedgraph_filtered[['chr', 'start', 'end', prefix]]

        dataframes.append(bedgraph_filtered) # Append DataFrames to List

    merged_dataframe = merge_dfs(dataframes)

    # Reorder columns to put control columns first
    control_cols = [col for col in merged_dataframe.columns if control in col]
    first_cols = ['chr', 'start', 'end'] + control_cols
    other_cols = [col for col in merged_dataframe.columns if col not in first_cols]
    merged_dataframe = merged_dataframe[first_cols + other_cols] 

    merged_dataframe = merged_dataframe.loc[:, ~merged_dataframe.columns.str.contains('total_reads')]

    for col in merged_dataframe.columns[3:]:
        merged_dataframe[['chr', 'start', 'end', col]].to_csv(f"{col}_library_trimmed_CpGs.bedGraph", sep='\t', index=False, header=False)

    return merged_dataframe

def comp_meth_frequencies(df, n_replicates, sample_list, control, cutoff, filter = False):
    """
    Filters regions based on cutoff for the difference between samples methylation frequency vs control methylation frequency.
    """
    
    included = pd.DataFrame() # DataFrame to store included regions
    filtered = pd.DataFrame() # DataFrame to store filtered regions

    for i in range(len(df)):
        count = 0

        for sample in sample_list:

            diff_samples = list() # List to store difference between samples
            diff_control = list() # List to store difference between controls

            for n in range(1, n_replicates + 1):
                
                if sample + "_rep" + str(n) in df.columns:
                    diff_samples.append(int(df[sample + "_rep" + str(n)].iloc[i])) # Append sample replicates to list
                else:
                    continue

                if control + "_rep" + str(n) in df.columns:
                    diff_control.append(int(df[control + "_rep" + str(n)].iloc[i])) # Append control replicates to list
                else:
                    continue   
            
            if abs(mean(diff_samples) - mean(diff_control)) >= cutoff:
                included = pd.concat([included, df.iloc[[i]]]) # Add differences in frequency >= cutoff between control and sample to included DataFrame 
                count += 1
            
            

        if count == 1 and filter == True:
            # If only one sample has difference in frequency >= cutoff between control and sample add to filtered DataFrame
            filtered = pd.concat([filtered, df.iloc[[i]]])
    
    included = included.drop_duplicates().reset_index(drop=True)
    filtered = filtered.drop_duplicates().reset_index(drop=True)

    return (included, filtered) if filter else included

def cutoff_heatmap(df, n_replicates, samples, n_samples, control, cutoff):
    """
    Filters regions for the hierarchical clustering.
    """
    samples_names = samples[-(int(n_samples) + 6):-6]

    included, filtered = comp_meth_frequencies(df, n_replicates, samples_names, control, cutoff, True)   

    return included, filtered

def cutoff_genm_regions(df, n_replicates, samples, n_samples, control, cutoff):
    """
    Filters regions for the genomic regions, for each sample.
    """
    samples_names = samples[-(int(n_samples) + 6):-6]

    for sample in samples_names:
        
        sample_df = comp_meth_frequencies(df, n_replicates, [sample], control, cutoff)
        sample_df.to_csv(sample +'__methylation.csv', sep='\t', index=False, header=False)

if __name__ == "__main__":
    # Merge bedGraph files and store in a file and trim with a cutoff of 5
    print(sys.argv)
    trim_bedgraph = trim_bedgraph(int(sys.argv[len(sys.argv)-1]), sys.argv[len(sys.argv)-5])

    trim_bedgraph.to_csv('correlation.csv', sep='\t', index=False)

    included_df, unfiltered_df = cutoff_heatmap(trim_bedgraph, int(sys.argv[len(sys.argv)-2]), list(sys.argv), sys.argv[len(sys.argv)-6], sys.argv[len(sys.argv)-5], int(sys.argv[len(sys.argv)-4]))

    removed_df, filtered_df = cutoff_heatmap(trim_bedgraph, int(sys.argv[len(sys.argv)-2]), list(sys.argv), sys.argv[len(sys.argv)-6], sys.argv[len(sys.argv)-5], int(sys.argv[len(sys.argv)-3]))

    # Save filtered and included regions to files
    filtered_df.to_csv('diff_methylation_filtered.csv', sep='\t', index=False, header=True)
    included_df.to_csv('diff_methylation_all.csv', sep='\t', index=False, header=True)

    cutoff_genm_regions(included_df, int(sys.argv[len(sys.argv)-2]), list(sys.argv), sys.argv[len(sys.argv)-6], sys.argv[len(sys.argv)-5], int(sys.argv[len(sys.argv)-4]))
