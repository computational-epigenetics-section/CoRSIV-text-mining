"""
Author: Wen-Jou Chang
Baylor College of Medicine

This script samples control regions from lookup tables based on CoRSIV metrics.
"""
import pandas as pd
import sys

# Get chromosome from command line argument
input_chr = sys.argv[1]

# Read CoRSIV metrics table and filter for specified chromosome
df = pd.read_csv('data/humanData/corsiv_control/corsiv_table_unique.csv') # unique combinations of CoRSIV metrics, can be obtained from process_control.ipynb
df = df[df["chr"]==input_chr]
grouped = df.groupby(by=['block_size'])

results = []

# Initialize empty dataframe to store control regions
df_list = pd.DataFrame(columns=["Region ID", "chr", "start", "end", "block_size", "CpG Count", "TSS Count", "GB Count", "TES Count", "EPIC Count", "HM450 Count", "Union Count", "gcount", "selected", "ID"])


# Function to filter out overlapping regions
def eliminate_overlapping_regions(df):
    if df.empty:
        return df

    # First pass: keep regions that don't overlap based on start/end coordinates
    keep_indices = []
    current_end = df.iloc[0]["end"]
    keep_indices.append(0)
    for i in range(1, len(df)):
        start, end = df.iloc[i]["start"], df.iloc[i]["end"]
        if start >= current_end:
            keep_indices.append(i)
            current_end = end
    df = df.iloc[keep_indices]
    
    # Second pass: check for overlaps with already selected regions in df_list
    def overlaps(row, ref_df):
        # Check if any rows in ref_df overlap with the given row
        overlap = ref_df[(ref_df['chr'] == row['chr']) & 
                            (ref_df['start'] < row['end']) & 
                            (ref_df['end'] > row['start'])]
        return not overlap.empty

    mask = df.apply(lambda row: not overlaps(row, df_list), axis=1)
    return df[mask]
    
# Process each block size group
for block_size_value, group in grouped:
    # Read control regions table for current block size
    temp = pd.read_csv(f"{input_chr}_{block_size_value}bp_control_table.csv")
    temp["gcount"] = temp['TSS Count'] +temp['GB Count'] + temp['TES Count']
    temp["selected"] = False

    # For each CoRSIV metric combination
    for index, row in group.iterrows():
        # Find matching control regions with exact metric matches
        matching_rows = temp[
        (temp["selected"] == False) &
        (temp['TSS Count'] == row['tss_count']) &
        (temp['GB Count'] == row['Gene_body_count']) &
        (temp['TES Count'] == row['tes_count']) &
        (temp['Union Count'] == row['Union Count']) &
        (temp['CpG Count'] == row['CpG_count'])
        ]
        iter_count = 0
        matching_rows = eliminate_overlapping_regions(matching_rows)

        # If not enough matches found, gradually relax matching criteria
        while matching_rows.shape[0] // row['count'] < 10:
            
            # First relaxation: Allow CpG count to vary by ±10%
            if iter_count == 0:
                matching_rows = temp[
                (temp["selected"] == False) &
                (temp['TSS Count'] == row['tss_count']) &
                (temp['GB Count'] == row['Gene_body_count']) &
                (temp['TES Count'] == row['tes_count']) &
                (temp['Union Count'] == row['Union Count']) &
                (int(row['CpG_count']*0.9) <= temp['CpG Count']) &
                (temp['CpG Count'] <= int(row['CpG_count']*1.1))
                ]
            # Second relaxation: Allow CpG count to vary by ±20%
            elif iter_count == 1:
                matching_rows = temp[
                (temp["selected"] == False) &
                (temp['TSS Count'] == row['tss_count']) &
                (temp['GB Count'] == row['Gene_body_count']) &
                (temp['TES Count'] == row['tes_count']) &
                (temp['Union Count'] == row['Union Count']) &
                (int(row['CpG_count']*0.8) <= temp['CpG Count']) &
                (temp['CpG Count'] <= int(row['CpG_count']*1.2))
                ]
            # Third relaxation: Match total gene feature count instead of individual counts
            elif iter_count == 2:
                ref_sum = row['tss_count']+row['Gene_body_count']+row['tes_count']
                matching_rows = temp[
                (temp["selected"] == False) &
                (ref_sum == temp["gcount"]) &
                (temp['Union Count'] == row['Union Count']) &
                (int(row['CpG_count']*0.8) <= temp['CpG Count']) &
                (temp['CpG Count'] <= int(row['CpG_count']*1.2))
                ]
            # Further relaxations: Allow gene feature count to vary by increasing percentages
            elif 3 <= iter_count <= 12:
                ref_sum = row['tss_count']+row['Gene_body_count']+row['tes_count']
                if ref_sum == 0:
                    matching_rows = temp[
                    (temp["selected"] == False) &
                    (0 <= temp["gcount"]) &
                    (temp["gcount"] <= 2) &
                    (temp['Union Count'] == row['Union Count']) &
                    (int(row['CpG_count']*0.8) <= temp['CpG Count']) &
                    (temp['CpG Count'] <= int(row['CpG_count']*1.2))
                    ]
                else:
                    matching_rows = temp[
                    (temp["selected"] == False) &
                    (int(ref_sum*(1-(iter_count-2)*0.1)) <= temp["gcount"]) &
                    (temp["gcount"] <= int(ref_sum*(1+(iter_count-2)*0.1))) &
                    (temp['Union Count'] == row['Union Count']) &
                    (int(row['CpG_count']*0.8) <= temp['CpG Count']) &
                    (temp['CpG Count'] <= int(row['CpG_count']*1.2))
                    ]
            # Final relaxations: Allow both gene features and CpG count to vary widely
            elif 13 <= iter_count <= 20:
                ref_sum = row['tss_count']+row['Gene_body_count']+row['tes_count']
                if ref_sum == 0:
                    matching_rows = temp[
                    (0 <= temp["gcount"]) &
                    (temp["gcount"] <= 2) &
                    (temp['Union Count'] == row['Union Count']) &
                    (int(row['CpG_count']*(2-iter_count*0.1)) <= temp['CpG Count']) &
                    (temp['CpG Count'] <= int(row['CpG_count']*(iter_count*0.1)))
                    ]
                else:
                    matching_rows = temp[
                    (temp["selected"] == False) &
                    (0 <= temp["gcount"]) &
                    (temp["gcount"] <= 2*ref_sum) &
                    (temp['Union Count'] == row['Union Count']) &
                    (int(row['CpG_count']*(2-iter_count*0.1)) <= temp['CpG Count']) &
                    (temp['CpG Count'] <= int(row['CpG_count']*(iter_count*0.1)))
                    ]
            else:
                break
            matching_rows = eliminate_overlapping_regions(matching_rows)
            iter_count += 1

        # Record results for this CoRSIV metric combination
        record = {
            'chr': input_chr,
            'block_size': block_size_value,
            'matching_count': matching_rows.shape[0],
            'CpG_count': row['CpG_count'],
            'tss_count': row['tss_count'],
            'Gene_body_count': row['Gene_body_count'],
            'tes_count': row['tes_count'],
            'Union Count': row['Union Count'],
            'CoRSIV Count': row['count']
        }
        results.append(record)

        # Add selected control regions to final list
        matching_rows["ID"] = [(row['CpG_count'], row['tss_count'], row['Gene_body_count'], row['tes_count'], row['Union Count'])] * len(matching_rows)
        if matching_rows.shape[0] // row['count'] <= 10:
            sampled_df_idx = matching_rows.index
        else:
            sampled_df_idx = matching_rows.sample(n=10*row["count"], replace=False).index
        temp.loc[sampled_df_idx, 'selected'] = True
        df_list = pd.concat([df_list, matching_rows.loc[sampled_df_idx]], ignore_index=True)

# Save results to files
results_df = pd.DataFrame(results)
df_list.to_csv(f"{input_chr}_control_candidates.csv", index=False) # contains all 10 sets of control regions, all files for all chromosomes are concatenated into control_candidates.csv
