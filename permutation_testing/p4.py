import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
from tqdm import tqdm
import numpy as np
import random
from collections import Counter
import time
from math import log10
from scipy import stats
from statistics import mean
from collections import defaultdict
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("disease", help="ex: anthropometric", default=None)

args = parser.parse_args()

disease = args.disease

corsiv_probe_list = pd.read_csv("Data/corsiv_all_probes.txt", sep="\t", header=None)
corsiv_probe_list = set(corsiv_probe_list.iloc[:, 3])
control_probe_list = pd.read_csv("Data/control_probes_4.txt", sep="\t", header=None)
control_probe_list = set(control_probe_list.iloc[:, 3])
epic = pd.read_csv("Data/EPIC.hg38.txt", sep="\t", header=None)
epic_probe_list = set(epic.iloc[:, 3])
hm450 = pd.read_csv("Data/HM450.hg38.txt", sep="\t", header=None)
hm450_probe_list = set(hm450.iloc[:, 3])
illumina = epic_probe_list.union(hm450_probe_list)
enrichment_cutoff = len(corsiv_probe_list) / len(illumina) * 1.5


def read_in_probes(disease):
    df = pd.read_csv("Data/%s_all_probes.csv" % disease, header=0)
    probe_list = df["probeId"].to_list()
    probe_list = [
        s for s in probe_list if (str(s).startswith("cg") or str(s).startswith("ch."))
    ]
    c = Counter(probe_list)
    c = {
        key: count
        for key, count in c.items()
        if (key in epic_probe_list) or (key in hm450_probe_list)
    }  # key being probe id, value being the number of papers it's found in
    return c


def get_ratio_pval(
    count_dictionary, shuffled_corsiv_probe_list, shuffled_control_probes_list
):
    """
    Given a dictionary where key is the probe and val is the number of papers a given probe shows up in,
    return enrichment ratio and p value.
    """
    paper_threshold_count = []
    corsiv_count = []
    control_count = []
    i = 1
    max_probe_count = max(
        count_dictionary.values()
    )  # count_dictionary stores probeid:number of papers, this line finds the highest number of papers reporting a probe
    while i <= max_probe_count:
        dummy_dict = {
            key: count for key, count in count_dictionary.items() if count == i
        }
        logval = log10(len(dummy_dict)) if len(dummy_dict) > 0 else 0
        paper_threshold_count.append(
            (i, logval)
        )  # i is number of papers, logval is the log10(nuber of probes reported)
        i += 1
    probe_cutoff = max_probe_count
    for i in range(paper_threshold_count[-1][0], 0, -1):
        # this for loop finds the highest number of papers reporting at least 10 probes,
        # we only include these points in our enrichment calculation
        if paper_threshold_count[i - 1][1] < 1:
            continue
        probe_cutoff = i
        break
    paper_threshold_count = paper_threshold_count[:probe_cutoff]
    i = 1
    while i <= probe_cutoff:
        dummy_dict = {
            key: count for key, count in count_dictionary.items() if count == i
        }
        corsiv_overlap_count = len(
            set(dummy_dict.keys()).intersection(shuffled_corsiv_probe_list)
        )
        corsiv_overlap = (
            corsiv_overlap_count / len(dummy_dict) * 100 if len(dummy_dict) > 0 else 0
        )
        corsiv_count.append((i, corsiv_overlap))
        control_overlap_count = len(
            set(dummy_dict.keys()).intersection(shuffled_control_probes_list)
        )
        control_overlap = (
            control_overlap_count / len(dummy_dict) * 100 if len(dummy_dict) > 0 else 0
        )
        control_count.append((i, control_overlap))
        i += 1
    corsiv_ratio = sum([i * pct for i, pct in corsiv_count])
    control_ratio = sum([i * pct for i, pct in control_count])
    if len(corsiv_count) == 1:
        enrichment_ratio = -1
    elif corsiv_ratio == 0 and control_ratio == 0:
        enrichment_ratio = -3
    elif control_ratio == 0:
        enrichment_ratio = -2
    else:
        enrichment_ratio = round(corsiv_ratio / control_ratio, 2)
    return enrichment_ratio


PATH = "%s_permutation_testing_results_chunks" % disease
if not os.path.exists(PATH):
    os.makedirs(PATH)

disease_probes_dict = read_in_probes(disease)


# 30001-31000


permutation_id_list = []
enrichment_ratio_list = []
p_value_list = []
disease_list = []

for i in range(30001, 31001):
    shuffled_df = pd.read_csv("permuted_sets/p_%s" % i, sep="\t", header=0)
    shuffled_corsiv_probe_list = list(
        set(list(shuffled_df.loc[shuffled_df["class"] == "corsiv"]["probeId"]))
    )
    shuffled_control_probe_list = list(
        set(list(shuffled_df.loc[shuffled_df["class"] == "control"]["probeId"]))
    )
    enrichment_ratio = get_ratio_pval(
        disease_probes_dict, shuffled_corsiv_probe_list, shuffled_control_probe_list
    )

    permutation_id_list.append("p_%s" % i)
    enrichment_ratio_list.append(enrichment_ratio)
    disease_list.append(disease)

final_df = pd.DataFrame(
    {
        "shuffled_set_id": permutation_id_list,
        "enrichment_ratio": enrichment_ratio_list,
        "Disease": disease_list,
    }
)

final_df.to_csv(
    "%s/%s_enrichment_after_permutations_%s.bed" % (PATH, disease, i),
    sep="\t",
    header=True,
    index=False,
)


# 31001-32000


permutation_id_list = []
enrichment_ratio_list = []
p_value_list = []
disease_list = []

for i in range(31001, 32001):
    shuffled_df = pd.read_csv("permuted_sets/p_%s" % i, sep="\t", header=0)
    shuffled_corsiv_probe_list = list(
        set(list(shuffled_df.loc[shuffled_df["class"] == "corsiv"]["probeId"]))
    )
    shuffled_control_probe_list = list(
        set(list(shuffled_df.loc[shuffled_df["class"] == "control"]["probeId"]))
    )
    enrichment_ratio = get_ratio_pval(
        disease_probes_dict, shuffled_corsiv_probe_list, shuffled_control_probe_list
    )

    permutation_id_list.append("p_%s" % i)
    enrichment_ratio_list.append(enrichment_ratio)
    disease_list.append(disease)

final_df = pd.DataFrame(
    {
        "shuffled_set_id": permutation_id_list,
        "enrichment_ratio": enrichment_ratio_list,
        "Disease": disease_list,
    }
)

final_df.to_csv(
    "%s/%s_enrichment_after_permutations_%s.bed" % (PATH, disease, i),
    sep="\t",
    header=True,
    index=False,
)


# 32001-33000


permutation_id_list = []
enrichment_ratio_list = []
p_value_list = []
disease_list = []

for i in range(32001, 33001):
    shuffled_df = pd.read_csv("permuted_sets/p_%s" % i, sep="\t", header=0)
    shuffled_corsiv_probe_list = list(
        set(list(shuffled_df.loc[shuffled_df["class"] == "corsiv"]["probeId"]))
    )
    shuffled_control_probe_list = list(
        set(list(shuffled_df.loc[shuffled_df["class"] == "control"]["probeId"]))
    )
    enrichment_ratio = get_ratio_pval(
        disease_probes_dict, shuffled_corsiv_probe_list, shuffled_control_probe_list
    )

    permutation_id_list.append("p_%s" % i)
    enrichment_ratio_list.append(enrichment_ratio)
    disease_list.append(disease)

final_df = pd.DataFrame(
    {
        "shuffled_set_id": permutation_id_list,
        "enrichment_ratio": enrichment_ratio_list,
        "Disease": disease_list,
    }
)

final_df.to_csv(
    "%s/%s_enrichment_after_permutations_%s.bed" % (PATH, disease, i),
    sep="\t",
    header=True,
    index=False,
)


# 33001-34000


permutation_id_list = []
enrichment_ratio_list = []
p_value_list = []
disease_list = []

for i in range(33001, 34001):
    shuffled_df = pd.read_csv("permuted_sets/p_%s" % i, sep="\t", header=0)
    shuffled_corsiv_probe_list = list(
        set(list(shuffled_df.loc[shuffled_df["class"] == "corsiv"]["probeId"]))
    )
    shuffled_control_probe_list = list(
        set(list(shuffled_df.loc[shuffled_df["class"] == "control"]["probeId"]))
    )
    enrichment_ratio = get_ratio_pval(
        disease_probes_dict, shuffled_corsiv_probe_list, shuffled_control_probe_list
    )

    permutation_id_list.append("p_%s" % i)
    enrichment_ratio_list.append(enrichment_ratio)
    disease_list.append(disease)

final_df = pd.DataFrame(
    {
        "shuffled_set_id": permutation_id_list,
        "enrichment_ratio": enrichment_ratio_list,
        "Disease": disease_list,
    }
)

final_df.to_csv(
    "%s/%s_enrichment_after_permutations_%s.bed" % (PATH, disease, i),
    sep="\t",
    header=True,
    index=False,
)


# 34001-35000


permutation_id_list = []
enrichment_ratio_list = []
p_value_list = []
disease_list = []

for i in range(34001, 35001):
    shuffled_df = pd.read_csv("permuted_sets/p_%s" % i, sep="\t", header=0)
    shuffled_corsiv_probe_list = list(
        set(list(shuffled_df.loc[shuffled_df["class"] == "corsiv"]["probeId"]))
    )
    shuffled_control_probe_list = list(
        set(list(shuffled_df.loc[shuffled_df["class"] == "control"]["probeId"]))
    )
    enrichment_ratio = get_ratio_pval(
        disease_probes_dict, shuffled_corsiv_probe_list, shuffled_control_probe_list
    )

    permutation_id_list.append("p_%s" % i)
    enrichment_ratio_list.append(enrichment_ratio)
    disease_list.append(disease)

final_df = pd.DataFrame(
    {
        "shuffled_set_id": permutation_id_list,
        "enrichment_ratio": enrichment_ratio_list,
        "Disease": disease_list,
    }
)

final_df.to_csv(
    "%s/%s_enrichment_after_permutations_%s.bed" % (PATH, disease, i),
    sep="\t",
    header=True,
    index=False,
)


# 35001-36000


permutation_id_list = []
enrichment_ratio_list = []
p_value_list = []
disease_list = []

for i in range(35001, 36001):
    shuffled_df = pd.read_csv("permuted_sets/p_%s" % i, sep="\t", header=0)
    shuffled_corsiv_probe_list = list(
        set(list(shuffled_df.loc[shuffled_df["class"] == "corsiv"]["probeId"]))
    )
    shuffled_control_probe_list = list(
        set(list(shuffled_df.loc[shuffled_df["class"] == "control"]["probeId"]))
    )
    enrichment_ratio = get_ratio_pval(
        disease_probes_dict, shuffled_corsiv_probe_list, shuffled_control_probe_list
    )

    permutation_id_list.append("p_%s" % i)
    enrichment_ratio_list.append(enrichment_ratio)
    disease_list.append(disease)

final_df = pd.DataFrame(
    {
        "shuffled_set_id": permutation_id_list,
        "enrichment_ratio": enrichment_ratio_list,
        "Disease": disease_list,
    }
)

final_df.to_csv(
    "%s/%s_enrichment_after_permutations_%s.bed" % (PATH, disease, i),
    sep="\t",
    header=True,
    index=False,
)


# 36001-37000


permutation_id_list = []
enrichment_ratio_list = []
p_value_list = []
disease_list = []

for i in range(36001, 37001):
    shuffled_df = pd.read_csv("permuted_sets/p_%s" % i, sep="\t", header=0)
    shuffled_corsiv_probe_list = list(
        set(list(shuffled_df.loc[shuffled_df["class"] == "corsiv"]["probeId"]))
    )
    shuffled_control_probe_list = list(
        set(list(shuffled_df.loc[shuffled_df["class"] == "control"]["probeId"]))
    )
    enrichment_ratio = get_ratio_pval(
        disease_probes_dict, shuffled_corsiv_probe_list, shuffled_control_probe_list
    )

    permutation_id_list.append("p_%s" % i)
    enrichment_ratio_list.append(enrichment_ratio)
    disease_list.append(disease)

final_df = pd.DataFrame(
    {
        "shuffled_set_id": permutation_id_list,
        "enrichment_ratio": enrichment_ratio_list,
        "Disease": disease_list,
    }
)

final_df.to_csv(
    "%s/%s_enrichment_after_permutations_%s.bed" % (PATH, disease, i),
    sep="\t",
    header=True,
    index=False,
)


# 37001-38000


permutation_id_list = []
enrichment_ratio_list = []
p_value_list = []
disease_list = []

for i in range(37001, 38001):
    shuffled_df = pd.read_csv("permuted_sets/p_%s" % i, sep="\t", header=0)
    shuffled_corsiv_probe_list = list(
        set(list(shuffled_df.loc[shuffled_df["class"] == "corsiv"]["probeId"]))
    )
    shuffled_control_probe_list = list(
        set(list(shuffled_df.loc[shuffled_df["class"] == "control"]["probeId"]))
    )
    enrichment_ratio = get_ratio_pval(
        disease_probes_dict, shuffled_corsiv_probe_list, shuffled_control_probe_list
    )

    permutation_id_list.append("p_%s" % i)
    enrichment_ratio_list.append(enrichment_ratio)
    disease_list.append(disease)

final_df = pd.DataFrame(
    {
        "shuffled_set_id": permutation_id_list,
        "enrichment_ratio": enrichment_ratio_list,
        "Disease": disease_list,
    }
)

final_df.to_csv(
    "%s/%s_enrichment_after_permutations_%s.bed" % (PATH, disease, i),
    sep="\t",
    header=True,
    index=False,
)


# 38001-39000


permutation_id_list = []
enrichment_ratio_list = []
p_value_list = []
disease_list = []

for i in range(38001, 39001):
    shuffled_df = pd.read_csv("permuted_sets/p_%s" % i, sep="\t", header=0)
    shuffled_corsiv_probe_list = list(
        set(list(shuffled_df.loc[shuffled_df["class"] == "corsiv"]["probeId"]))
    )
    shuffled_control_probe_list = list(
        set(list(shuffled_df.loc[shuffled_df["class"] == "control"]["probeId"]))
    )
    enrichment_ratio = get_ratio_pval(
        disease_probes_dict, shuffled_corsiv_probe_list, shuffled_control_probe_list
    )

    permutation_id_list.append("p_%s" % i)
    enrichment_ratio_list.append(enrichment_ratio)
    disease_list.append(disease)

final_df = pd.DataFrame(
    {
        "shuffled_set_id": permutation_id_list,
        "enrichment_ratio": enrichment_ratio_list,
        "Disease": disease_list,
    }
)

final_df.to_csv(
    "%s/%s_enrichment_after_permutations_%s.bed" % (PATH, disease, i),
    sep="\t",
    header=True,
    index=False,
)


# 39001-40000


permutation_id_list = []
enrichment_ratio_list = []
p_value_list = []
disease_list = []

for i in range(39001, 40001):
    shuffled_df = pd.read_csv("permuted_sets/p_%s" % i, sep="\t", header=0)
    shuffled_corsiv_probe_list = list(
        set(list(shuffled_df.loc[shuffled_df["class"] == "corsiv"]["probeId"]))
    )
    shuffled_control_probe_list = list(
        set(list(shuffled_df.loc[shuffled_df["class"] == "control"]["probeId"]))
    )
    enrichment_ratio = get_ratio_pval(
        disease_probes_dict, shuffled_corsiv_probe_list, shuffled_control_probe_list
    )

    permutation_id_list.append("p_%s" % i)
    enrichment_ratio_list.append(enrichment_ratio)
    disease_list.append(disease)

final_df = pd.DataFrame(
    {
        "shuffled_set_id": permutation_id_list,
        "enrichment_ratio": enrichment_ratio_list,
        "Disease": disease_list,
    }
)


final_df.to_csv(
    "%s/%s_enrichment_after_permutations_%s.bed" % (PATH, disease, i),
    sep="\t",
    header=True,
    index=False,
)

