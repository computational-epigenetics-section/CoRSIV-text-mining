"""
Author: Wen-Jou Chang
Baylor College of Medicine

This script generates one lookup table for control regions per given region size on a given chromosome. 
Control regions are later randomly sampled from the respective lookup table based on CoRSIV metrics.
"""

import pandas as pd
import numpy as np
import sys

# parameters
chr = sys.argv[1] 

# bed file for gene bodies, promoters (padded 3kb), and three prime regions (padded 3kb), all obtained from UCSC
gb = pd.read_csv('hg38_gene_bodies.tab.txt', header=None, names=['chr', 'start', 'end', 'ensg', 'gene'], sep="\t")
gb = gb[gb['chr'] == chr]

fivep = pd.read_csv('hg38_promoters_3kb.tab.txt', header=None, names=['chr', 'start', 'end', 'ensg', 'gene'], sep="\t")
fivep = fivep[fivep['chr'] == chr]

threep = pd.read_csv('hg38_three_prime_region_3kb.tab.txt', header=None, names=['chr', 'start', 'end', 'ensg', 'gene'], sep="\t")
threep = threep[threep['chr'] == chr]

# 100bp bin-level annotation of CpG counts, TSS counts, GB counts, EPIC probe counts, bins overlapping CoRSIVs removed 
bins = pd.read_csv(f"{chr}_no_corsiv.csv")

# bed file for EPIC probes
epic = pd.read_csv('EPIC.hg38.txt', header=None, names=['chr', 'start', 'end', 'id'], sep="\t")
epic = epic[epic['chr'] == chr]

# bed file for HM450 probes
hm450 = pd.read_csv('HM450.hg38.txt', header=None, names=['chr', 'start', 'end', 'id'], sep="\t")
hm450 = hm450[hm450['chr'] == chr]

# bed file for EPIC_HM450 unified probe counts
probe = pd.read_csv('EPIC_HM450_control_table.bed', header=None, names=['chr', 'start', 'end', 'id', 'pcount'], sep="\t")
probe = probe[probe['chr'] == chr]

# bed file for annotated CoRSIV regions
corsiv = pd.read_csv("annotated_corsiv_all.csv")
bp_set = set(corsiv[corsiv["chr"] == chr]["block_size"])

def get_gene_count(region):
    gb_count = ((gb['start'] < region['end']) & (gb['end'] > region['start'])).sum()
    tss_count = ((fivep['start'] < region['end']) & (fivep['end'] > region['start'])).sum()
    tes_count = ((threep['start'] < region['end']) & (threep['end'] > region['start'])).sum()
    epic_count = ((epic['start'] < region['end']) & (epic['end'] > region['start'])).sum()
    hm450_count = ((hm450['start'] < region['end']) & (hm450['end'] > region['start'])).sum()
    probe_count = probe[(probe['start'] < region['end']) & (probe['end'] > region['start'])]['pcount'].sum()
    return tss_count, gb_count, tes_count, epic_count, hm450_count, probe_count

def build_lookup(block_size):
    n = block_size // 100
    records = []
    for i in range(0, len(bins) - n + 1):
        curr = bins.iloc[i:i + n, :]
        s, e = curr.iloc[0]['start'], curr.iloc[-1]['end']
        cpg_ct = curr['CpG Count'].sum()
        if e - s > n * 100:
            continue
        region = {'start': s, 'end': e}
        tss_ct, gb_ct, tes_ct, epic_ct, hm450_ct, probe_ct = get_gene_count(region)
        record = {
            'Region ID': f'{chr}_{s}_{e}',
            'chr': chr,
            'start': s,
            'end': e,
            'block_size': e - s,
            'CpG Count': cpg_ct,
            'TSS Count': tss_ct,
            'GB Count': gb_ct,
            'TES Count': tes_ct,
            'EPIC Count': epic_ct,
            'HM450 Count': hm450_ct,
            'Union Count': probe_ct
        }
        records.append(record)

    df = pd.DataFrame(records)
    df.to_csv(f'{chr}_{block_size}bp_control_table.csv', index=False)

for bp in bp_set:
    build_lookup(bp)
