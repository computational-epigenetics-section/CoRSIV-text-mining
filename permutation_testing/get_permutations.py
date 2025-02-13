import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
from tqdm import tqdm
import numpy as np
import random

PATH = "permuted_sets"
if not os.path.exists(PATH):
    os.makedirs(PATH)

for i in range(10):
    corsiv_all_probes = pd.read_csv("Data/corsiv_all_probes.txt", sep="\t", header=None)
    control_all_probes = pd.read_csv(
        "Data/control_probes_%s.txt" % (i + 1), sep="\t", header=None
    )

    corsiv_all_probes_modified = corsiv_all_probes.copy()
    control_all_probes_modified = control_all_probes.copy()
    control_all_probes_modified = control_all_probes_modified[[0, 1, 2, 3]]
    corsiv_all_probes_modified["class"] = "corsiv"
    control_all_probes_modified["class"] = "control"
    corsiv_all_probes_modified.columns = ["chr", "start", "end", "probeId", "class"]
    control_all_probes_modified.columns = ["chr", "start", "end", "probeId", "class"]

    concat_probes = pd.concat([corsiv_all_probes_modified, control_all_probes_modified])
    concat_probes = concat_probes.reset_index(drop=True)
    concat_probes = concat_probes[["probeId", "class"]]

    initial_order = list(concat_probes["class"])

    order_list = []
    order_list.append(initial_order)

    # permutation_count = 100001
    count = 0
    for j in range((i * 10000), ((i + 1) * 10000)):
        new_order = random.sample(initial_order, len(initial_order))
        if new_order in order_list:
            continue
        else:
            order_list.append(new_order)
            temp_df = pd.DataFrame(
                {"probeId": list(concat_probes["probeId"]), "class": new_order}
            )
        temp_df.to_csv("%s/p_%s" % (PATH, (j + 1)), sep="\t", header=True, index=False)
        count += 1

    print("Permutations with control set %s : %s", (i, count))
