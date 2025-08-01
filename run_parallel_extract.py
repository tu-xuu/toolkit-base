"""
Parallel extractor for pollutant and meteorological NetCDF data.

This script uses multiprocessing to call `extract_loc_all_data()` 
and extract daily pollutant and meteorological data (2023–2024) 
for each specified location (from loc_info.csv), 
while skipping IDs already processed.

Paths and parameters should be modified according to your setup.

Author: Chengxu Tong
"""


import os
import pandas as pd
from multiprocessing import Pool
from Code import extract_loc_all_data

def run_extract_for_id(loc_id):
    extract_loc_all_data(
        loc_id=loc_id,
        pollutants='all',
        pollutant_years=['2023','2024'],
        meteorology_years=['2023', '2024'],
        loc_csv_path="/path/to/your/city_data/loc_info.csv",
        output_path="/path/to/your/city_data/Results"
    )

def get_all_loc_ids(loc_csv_path):
    loc_df = pd.read_csv(os.path.expanduser(loc_csv_path))
    return loc_df['loc_id'].unique().tolist()

def get_done_ids(results_dir):
    results_dir = os.path.expanduser(results_dir)
    done_ids = []
    for fname in os.listdir(results_dir):
        if fname.startswith("loc_id_") and fname.endswith(".csv"):
            try:
                loc_id = int(fname.split("_")[2].split(".")[0])
                done_ids.append(loc_id)
            except:
                continue
    return sorted(set(done_ids))

if __name__ == "__main__":
    all_ids = get_all_loc_ids("/path/to/your/city_data/loc_info.csv")
    done_ids = get_done_ids("/path/to/your/city_data/Results")
    remaining_ids = [i for i in all_ids if i not in done_ids]

    print(f"All: {len(all_ids)} IDs")
    print(f"Done: {len(done_ids)} IDs")
    print(f"Remaining: {len(remaining_ids)} IDs")

    with Pool(processes=90) as pool:
        pool.map(run_extract_for_id, remaining_ids)