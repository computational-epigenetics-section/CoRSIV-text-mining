import sys
import pandas as pd

chr = sys.argv[1]
bp = sys.argv[2]
gcount = int(sys.argv[3])
pcount = int(sys.argv[4])
cpg = int(sys.argv[5])
cpg_ratio = float(sys.argv[6])

def eliminate_overlapping_regions(df):
    if df.empty:
        return df
    # df.sort_values(by='start', inplace=True)
    keep_indices = []
    current_end = df.iloc[0]["end"]
    keep_indices.append(0)
    for i in range(1, len(df)):
        start, end = df.iloc[i]["start"], df.iloc[i]["end"]
        if start >= current_end:
            keep_indices.append(i)
            current_end = end
    df = df.iloc[keep_indices].reset_index(drop=True)
    return df

temp = pd.read_csv(f"/storage/waterland/home/u239646/text_mining/control/table/chr{chr}_{bp}bp_control_table.csv")
temp["gcount"] = temp['TSS Count'] +temp['GB Count'] + temp['TES Count']
matching_rows = temp[(0 <= temp["gcount"]) & (temp["gcount"] <= gcount*2) & (temp['Union Count'] == pcount) & (int(cpg*(1.0-cpg_ratio)) <= temp['CpG Count']) & (temp['CpG Count'] <= int(cpg*(1.0+cpg_ratio)))]
matching_rows = eliminate_overlapping_regions(matching_rows)
print(matching_rows)