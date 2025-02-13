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

PATH = "final_permutation_testing_results"
if not os.path.exists(PATH):
    os.makedirs(PATH)

overall_results = pd.DataFrame()
for i in range(1000, 100001, 1000):
    temp = pd.read_csv(
        "%s_permutation_testing_results_chunks/%s_enrichment_after_permutations_%s.bed"
        % (disease, disease, i),
        sep="\t",
        header=0,
    )
    overall_results = pd.concat([overall_results, temp])


print(overall_results.shape[0])

overall_results["shuffled_set_id_int"] = overall_results["shuffled_set_id"].str.split("_").str[-1]

overall_results["shuffled_set_id_int"] = overall_results["shuffled_set_id_int"].astype(
    int
)
overall_results = overall_results.sort_values("shuffled_set_id_int")
overall_results = overall_results[["shuffled_set_id", "enrichment_ratio", "Disease"]]
overall_results.to_csv(
    "%s/%s_enrichment_after_permutations.bed"
    % (PATH, disease),
    sep="\t",
    header=True,
    index=False,
)
