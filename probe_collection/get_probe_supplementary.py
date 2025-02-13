import pandas as pd
import re
from collections import defaultdict
import os
import sys
import pandas.io.common

probe_pattern = r'cg\d{8}'
probe_pattern2 = r'ch\.(?:[1-9]|1[0-9]|2[0-2])\.\d+[A-Z]'
input_cat = sys.argv[1]

def get_csv_probes(format):
    folder = f"YOUR_DIRECTORY/{input_cat}/{format}"
    filenames = os.listdir(folder)
    filenames = [f for f in filenames if f.endswith("csv")]
    invalid_name = []
    cg_probes = defaultdict(set)
    for name in filenames:
        study_id = name.split("-")[0]
        try:
            df = pd.read_csv(f"{folder}/{name}")
        except UnicodeDecodeError:
            df = pd.read_csv(f"{folder}/{name}", encoding='cp1252')
        except pandas.errors.EmptyDataError:
            print(f"{name} is empty.")
            continue
        except:
            try:
                df = pd.read_csv(f"{folder}/{name}", sep=None, engine="python", error_bad_lines=False)
            except:
                invalid_name.append(name)
                continue
        if df.empty:
            continue

        df = df.to_string()
        probe = re.findall(probe_pattern, df)
        probe2 = re.findall(probe_pattern2, df)
        if len(probe) + len(probe2) <= 1000:
            cg_probes[study_id].update(probe)
            cg_probes[study_id].update(probe2)
    data = [{'pmcid': key, 'probeId': probe_id} for key, value in cg_probes.items() for probe_id in value]
    df = pd.DataFrame(data)
    df.to_csv(f"probe/{input_cat}_{format}_probes.csv", index=0)
    print(f"[{input_cat} ({format})] Skipped files: {invalid_name}")

def get_excel_probes(format):
    folder = f"{input_cat}/{format}"
    filenames = os.listdir(folder)
    filenames = [f for f in filenames if f.endswith(format)]
    invalid_name = []
    cg_probes = defaultdict(set)
    for k in range(len(filenames)):
        name = filenames[k]
        study_id = name.split("-")[0]
        # print(k, name)
        try:
            xl = pd.ExcelFile(f"{folder}/{name}")
        except:
            invalid_name.append(name)
            continue
        sheet_names = xl.sheet_names  # see all sheet names
        for n in sheet_names:
            try:
                df = xl.parse(n)  # read a specific sheet to DataFrame
            except:
                invalid_name.append(f"{name}-{sheet_names}")
                continue
            if df.empty:
                continue
            df = df.to_string()
            probe = re.findall(probe_pattern, df)
            probe2 = re.findall(probe_pattern2, df)
            if len(probe) + len(probe2) <= 1000: # limit the number of probes to 1000
                cg_probes[study_id].update(probe)
                cg_probes[study_id].update(probe2)
    data = [{'pmcid': key, 'probeId': probe_id} for key, value in cg_probes.items() for probe_id in value]
    df = pd.DataFrame(data)
    df.to_csv(f"YOUR_DIRECTORY/probe/{input_cat}_{format}_probes.csv", index=0)
    print(f"[{input_cat}] Skipped files: {invalid_name}")

get_csv_probes("csv")
get_csv_probes("xls")
get_csv_probes("xlsx")
dirs = os.listdir(f"YOUR_DIRECTORY/text_mining/{input_cat}")
excel_formats = [d for d in dirs if (d.lower()== "xlsm" or d.lower()== "xlsb")]
for l in excel_formats:
    get_excel_probes(l)