import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from math import log10
from statistics import mean
import os
from matplotlib import ticker


PROJECT_PATH = "YOUR DIRECTORY"
if PROJECT_PATH:
    os.chdir(PROJECT_PATH)
CATEGORIES = [["Neoplasms"], ["Cardiovascular Diseases"], ["Digestive System Diseases"], ["Endocrine System Diseases"], ["Hemic and Lymphatic Diseases"], ["Immune System Diseases"], ["Metabolic Diseases"], ["Mental Disorders", "Nervous System Diseases"], ["Obesity"], ["Respiratory Tract Diseases"], ["Urogenital Diseases"]]
CATEGORY_NAMES = ["cancer", "cardiovascular", "digestive", "endocrine", "hematological", "immune", "metabolic", "neurological", "obesity", "respiratory", "urogenital"]#, "T2D"]
COLOR_TEMPLATE = ['#e6194B','#f58231','#f3c300','#469990','#808000','#2f8e3b','#0db7dd','#4363d8','#800000','#f032e6','#911eb4']#, "#000075"]#'#8298e5'
CONTROLS = []
CORSIV_PROBE_DF = pd.read_csv("data/humanData/corsiv_control/corsiv_all_probes_id.txt", sep="\t", names=["chr", "start", "end", "probeId", "corsiv_start", "corsiv_end", "corsiv_id"])
CORSIV_PROBE_LIST = set(CORSIV_PROBE_DF.iloc[:,3])
for i in range(1, 11):
    control_probe_df = pd.read_csv(f"data/humanData/corsiv_control/control_probes_{i}.txt", sep="\t", names=["chr", "start", "end", "probeId", "_", "corsiv_start", "corsiv_end", "corsiv_id"])
    control_probe_list = set(control_probe_df.iloc[:,3])
    CONTROLS.append(control_probe_list)
    
def read_in_probes(input_cat):
    cat = "metabolic_diseases" if input_cat == "metabolic" else input_cat
    df = pd.read_csv(f"data/probe/{cat}_all_probes.csv")
    probe_list = df["probeId"].to_list()
    c = dict(Counter(probe_list))
    return c


def calculate_points(count_dictionary, input_cat):
    paper_threshold_count = []
    corsiv_count = []
    control_count = []
    i = 1
    if not count_dictionary:
        return [], [], [], 0, 0, set([])
    max_probe_count = max(count_dictionary.values())
    probe_res_count = 0
    corsiv_paper_set = set()
    while i <= max_probe_count:
        dummy_dict = {key:count for key, count in count_dictionary.items() if count == i}
        logval = len(dummy_dict)
        paper_threshold_count.append((i, logval))
        i += 1
    probe_cutoff = max_probe_count
    for i in range(paper_threshold_count[-1][0], 0, -1):
        if paper_threshold_count[i-1][1] < 10:
            continue
        probe_cutoff = i
        break
    paper_threshold_count = paper_threshold_count[:probe_cutoff]
    i = 1
    if isinstance(input_cat, str):
        cat = "metabolic_diseases" if input_cat == "metabolic" else input_cat
        df = pd.read_csv(f"probe/{cat}_all_probes.csv")
    else:
        df = input_cat
    while i <= probe_cutoff:
        dummy_dict = {key:count for key, count in count_dictionary.items() if count == i}
        overlapping_probes = set(dummy_dict.keys()).intersection(CORSIV_PROBE_LIST)
        corsiv_overlap_count = len(overlapping_probes)
        corsiv_overlap = corsiv_overlap_count / len(dummy_dict) *100 if len(dummy_dict) > 0 else 0
        corsiv_count.append((i, corsiv_overlap))
        curr_control = []
        for control_set in CONTROLS:
            control_overlap_count = len(set(dummy_dict.keys()).intersection(control_set))
            control_overlap = control_overlap_count / len(dummy_dict) *100 if len(dummy_dict) > 0 else 0
            curr_control.append(control_overlap)
        control_count.append((i, mean(curr_control)))
        if corsiv_overlap > mean(curr_control): #and i > 1:
            probe_res_count += corsiv_overlap_count
            corsiv_paper_set.update(df[df["probeId"].isin(overlapping_probes)]["pmcid"].unique())
        i += 1
    return paper_threshold_count, corsiv_count, control_count, probe_res_count, len(corsiv_paper_set), corsiv_paper_set

def plot_enrichment(counts, title_text, cat_index, output=None, show_ratio=True, format="svg", show_legend=False, show_figure=True, show_y_label=True, box_placement=(0.15, 0.9), export_all=False):
    if counts[0] == []:
        return 0, 0
    fig = plt.figure(figsize=(6, 6))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 2])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)     
    
    x_values, y_values = zip(*counts[0])
    ax1.plot(x_values, y_values, marker='o', linestyle='-', color = COLOR_TEMPLATE[cat_index])
    ax1.set_title(title_text.capitalize() if " " not in title_text else title_text, fontsize=30)# .capitalize()
    max_y = int(log10(max(y_values))) + 1
    yticks = [10**i for i in range(1, max_y + 1)]
    ax1.set_yscale('log')
    ax1.set_yticks(yticks)
    ax1.set_yticklabels([f'$10^{i}$' for i in range(1, max_y + 1)])
    # ax1.yaxis.set_minor_locator(ticker.LogLocator(subs=range(2, 10)))
    ax1.set_xticks(range(1, int(max(x_values))+1))
    ax1.tick_params(axis='x', bottom=True, direction='inout', labelbottom=False, length=10)
    ax1.tick_params(axis='y', which='both', left=True, labelleft=True)
    x_values, y_values = zip(*counts[1])
    ax2.plot(x_values, y_values, marker='o', linestyle='-', color = COLOR_TEMPLATE[cat_index], label="CoRSIV")
    x_values, y_values = zip(*counts[2])
    ax2.plot(x_values, y_values, marker='o', linestyle='-', color = "grey", label="Control")
    ax2.set_xticks(range(1, int(max(x_values))+1))
    ax2.set_xlabel('Number of Papers Reporting Probe', fontsize=18)
    if show_y_label:
        ax1.set_ylabel('Number of Probes', fontsize=16)
        ax2.set_ylabel('Overlapping Probes (%)', fontsize=16)
        ax2.tick_params(axis='y', labelsize=18)
    ax2.tick_params(axis='x', which='both', bottom=True, labelbottom=True, labelsize=18)
    if not show_y_label:
        decimals = 1 if any(tick % 1 != 0 for tick in ax2.get_yticks()) else 0
        ax2.yaxis.set_major_formatter(ticker.PercentFormatter(decimals=decimals))
    enrichment_ratio = 0
    if show_ratio:
        corsiv_ratio = sum([i*pct for i, pct in counts[1]])
        control_ratio = sum([i*pct for i, pct in counts[2]])
        enrichment_ratio = round(corsiv_ratio / control_ratio, 1) if control_ratio != 0 else -1
        middle_text = f"Ratio = {enrichment_ratio}\n{counts[3]:,} Probes\n{counts[4]:,} Papers" if enrichment_ratio > 1 else "No Enrichment"
        
        ax2.annotate(middle_text,
                    xy=box_placement,  # Adjust position to account for bbox size
                    xycoords='axes fraction',  # Use axes fraction for coordinates
                    ha='center',  # Center the text horizontally
                    va='top',     # Align text to the top of the box
                    fontsize=16,  # Font size
                    color="red",
                    bbox=dict(facecolor='white', edgecolor='red', linewidth=2, boxstyle='square,pad=0.5'))  # Add thicker square box around text
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    
    if show_legend:
        ax2.legend(bbox_to_anchor=(0.35, 0.5), shadow=False, frameon=False)

    if output is None:
        if show_figure:
            plt.show() 
    else:
        if export_all or (enrichment_ratio > 1 and len(counts[0]) > 1):
            fig.savefig(output, format=format, bbox_inches='tight')
    plt.close()
    return len(counts[0]), enrichment_ratio


def breakdown(df, term_pmcid_map, disease, cat_index, output=None, show_figure=False, show_y_label=False, show_legend=False, format="svg", box_placement=(0.15, 0.9), export_all=False):
    df_categorized = df[df["pmcid"].isin(set().union(*[term_pmcid_map[d] for d in disease]))]
    probe_list = df_categorized["probeId"].to_list()
    c = Counter(probe_list)
    l1, l2, l3, p, p2, paper_set = calculate_points(c, df_categorized)
    paper, r = plot_enrichment([l1, l2, l3, p, p2], disease[0] if len(disease) == 1 else "Neurological", cat_index, output=output, format=format, show_figure=show_figure, show_legend=show_legend, show_y_label=show_y_label, box_placement=box_placement, export_all=export_all)
    return p, paper, r, p2, paper_set


def export_paper(broad_cat, diseases_name):
    cat = "metabolic_diseases" if broad_cat == "metabolic" else broad_cat
    neuro_df = pd.read_csv(f"probe/{cat}_all_probes.csv")
    if broad_cat == "metabolic":
        neuro_df["Filtered Mesh Term"] = neuro_df["Filtered Mesh Term"].apply(eval)
    else:
        neuro_df["Filtered Mesh Term"] = neuro_df["Filtered Mesh Term"].apply(lambda x: [term.strip() for term in x.split("|")])
    neuro_df = neuro_df[neuro_df["Filtered Mesh Term"].apply(lambda x: diseases_name in x) & (neuro_df["probeId"].isin(CORSIV_PROBE_LIST))]
    neuro_df = neuro_df.groupby("pmcid").agg({"probeId": lambda x: ";".join(x)}).reset_index()
    papers = pd.read_csv(f"pubmed_search/{cat}_final.csv")
    m = pd.merge(neuro_df, papers, left_on="pmcid", right_on="PMCID", how="left")
    m["Link to Paper"] = m["PMID"].apply(lambda x: f"https://pubmed.ncbi.nlm.nih.gov/{x}")
    m.rename(columns={"probeId": "CoRSIV Probes"}, inplace=True)
    m = m[["PMID", "Last Name", "Year", "Journal", "Title", "Abstract", "CoRSIV Probes", "Link to Paper"]]
    return m