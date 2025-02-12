import pandas as pd
import numpy as np
import sys

# parameters
chr = sys.argv[1] 

gb = pd.read_csv('/storage/waterland/home/u239646/text_mining/control/data/hg38_gene_bodies.tab.txt', header=None, names=['chr', 'start', 'end', 'ensg', 'gene'], sep="\t")
gb = gb[gb['chr'] == chr]

fivep = pd.read_csv('/storage/waterland/home/u239646/text_mining/control/data/hg38_promoters_3kb.tab.txt', header=None, names=['chr', 'start', 'end', 'ensg', 'gene'], sep="\t")
fivep = fivep[fivep['chr'] == chr]

threep = pd.read_csv('/storage/waterland/home/u239646/text_mining/control/data/hg38_three_prime_region_3kb.tab.txt', header=None, names=['chr', 'start', 'end', 'ensg', 'gene'], sep="\t")
threep = threep[threep['chr'] == chr]

bins = pd.read_csv(f"/storage/waterland/home/u239646/text_mining/control/data/{chr}_no_corsiv.csv")

epic = pd.read_csv('/storage/waterland/home/u239646/text_mining/control/data/EPIC.hg38.txt', header=None, names=['chr', 'start', 'end', 'id'], sep="\t")
epic = epic[epic['chr'] == chr]

hm450 = pd.read_csv('/storage/waterland/home/u239646/text_mining/control/data/HM450.hg38.txt', header=None, names=['chr', 'start', 'end', 'id'], sep="\t")
hm450 = hm450[hm450['chr'] == chr]

probe = pd.read_csv('/storage/waterland/home/u239646/text_mining/control/data/EPIC_HM450_control_table.bed', header=None, names=['chr', 'start', 'end', 'id', 'pcount'], sep="\t")
probe = probe[probe['chr'] == chr]

corsiv = pd.read_csv("/storage/waterland/home/u239646/text_mining/control/data/annotated_corsiv_all.csv")
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
    df.to_csv(f'/storage/waterland/home/u239646/text_mining/control/table/{chr}_{block_size}bp_control_table.csv', index=False)

for bp in bp_set:
    build_lookup(bp)
