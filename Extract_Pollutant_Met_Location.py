#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract daily pollutant (provided by Jing Wei) and meteorological data (ERA5) from NetCDF files for specified lat/lon locations, and merge into yearly CSVs.

Supports:
- Auto-detection of pollutants or manual specification
- Append all meteorological parameters from NC file
- Optional filtering by location ID
- Day-By-Day data batch processing
- Output as tidy CSVs for further use

Author: Chengxu Tong
"""

import os
import pandas as pd
import numpy as np
from netCDF4 import Dataset
from glob import glob
from tqdm import tqdm


def extract_pollutants(nc_root_path, loc_csv_path, pollutants, years, output_path, loc_id_list=None):
    """
    Extract specified pollutant(s) for selected loc points and save to CSV, year by year.

    Parameters:
        nc_root_path (str): Root folder containing pollutant-year subfolders
        loc_csv_path (str): Path to loc_info.csv
        pollutants (list or str): Pollutant variable(s) to extract, or 'all' to auto-detect
        years (list): Years to extract (e.g., ['2023', '2024'])
        output_path (str): Folder to save output CSVs
        loc_id_list (list or None): Optional list of loc IDs to extract; None = all
    """
    os.makedirs(output_path, exist_ok=True)
    nc_root_path = os.path.expanduser(nc_root_path)
    loc_csv_path = os.path.expanduser(loc_csv_path)
    output_path = os.path.expanduser(output_path)

   
    if pollutants == 'all':
        pollutants = []
        for item in os.listdir(nc_root_path):
            if os.path.isdir(os.path.join(nc_root_path, item)):
                parts = item.split("_")
                if len(parts) == 2:
                    pollutant_name = parts[0]
                    if pollutant_name not in pollutants:
                        pollutants.append(pollutant_name)
        print(f"[Auto] Detected pollutants: {pollutants}")

   
    loc_df = pd.read_csv(loc_csv_path)
    loc_df['longitude_wgs'] = loc_df['longitude_wgs'].astype(float)
    loc_df['latitude_wgs'] = loc_df['latitude_wgs'].astype(float)

 
    if loc_id_list is not None:
        loc_df = loc_df[loc_df['loc_id'].isin(loc_id_list)]

    for year in years:
        print(f"\n[Processing Year] {year}")

        loc_data = {
            loc_id: {
                'loc_id': loc_id,
                'longitude_wgs': row['longitude_wgs'],
                'latitude_wgs': row['latitude_wgs'],
                'records': []
            }
            for loc_id, row in loc_df.set_index('loc_id').iterrows()
        }

        for pollutant in pollutants:
            subfolder = f"{pollutant}_{year}"
            folder_path = os.path.join(nc_root_path, subfolder)
            if not os.path.exists(folder_path):
                print(f"[Skip] Folder not found: {subfolder}")
                continue

            nc_files = sorted(glob(os.path.join(folder_path, f"CHAP_{pollutant}_D1K_*_V*.nc")))
            if not nc_files:
                print(f"[Warning] No files found in {subfolder}")
                continue

            for file_path in tqdm(nc_files, desc=f"{pollutant}_{year}"):
                filename = os.path.basename(file_path)
                try:
                    date_str = filename.split("_")[3]
                    date = pd.to_datetime(date_str, format="%Y%m%d")
                except Exception as e:
                    print(f"[Error] Failed to parse date from {filename}: {e}")
                    continue

                try:
                    ds = Dataset(file_path)
                    data = np.array(ds[pollutant][:], dtype=np.float32)
                    fill_value = ds[pollutant]._FillValue
                    data[data == fill_value] = np.nan
                    lon = np.array(ds['lon'][:])
                    lat = np.array(ds['lat'][:])
                except Exception as e:
                    print(f"[Error] Failed to read {file_path}: {e}")
                    continue

                for loc_id in loc_data:
                    lon0 = loc_data[loc_id]['longitude_wgs']
                    lat0 = loc_data[loc_id]['latitude_wgs']
                    lon_idx = np.abs(lon - lon0).argmin()
                    lat_idx = np.abs(lat - lat0).argmin()
                    val = data[lat_idx, lon_idx]
                    val = None if np.isnan(val) else val

                    existing = next((r for r in loc_data[loc_id]['records'] if r['date'] == date), None)
                    if existing:
                        existing[pollutant] = val
                    else:
                        loc_data[loc_id]['records'].append({'date': date, pollutant: val})

                ds.close()

       
        for loc_id, info in loc_data.items():
            records = info['records']
            if not records:
                continue
            df = pd.DataFrame(records)
            df['loc_id'] = loc_id
            df['longitude_wgs'] = info['longitude_wgs']
            df['latitude_wgs'] = info['latitude_wgs']

            base_cols = ['loc_id', 'longitude_wgs', 'latitude_wgs', 'date']
            pollutant_cols = [col for col in pollutants if col in df.columns]
            df = df[base_cols + pollutant_cols]

            df.sort_values('date', inplace=True)
            output_file = os.path.join(output_path, f"loc_id_{loc_id}_{year}.csv")
            df.to_csv(output_file, index=False)
            print(f"[Done] Saved {os.path.basename(output_file)}")
            
def extract_meteorology(nc_file_paths, loc_csv_path, variables='all', output_path='.', loc_id_list=None):
    """
    Extract daily meteorological variables for selected loc points and append to existing pollutant CSVs.

    Parameters:
        nc_file_paths (dict): Mapping from year (str) to .nc file path, e.g., {'2023': 'ERA52023_daily.nc'}
        loc_csv_path (str): Path to loc_info.csv
        variables (list or 'all'): Meteorological variables to extract, or 'all' to auto-detect from .nc file
        output_path (str): Folder where loc_id_*.csv already exist and should be updated
        loc_id_list (list or None): If provided, only process these loc IDs
    """
    loc_csv_path = os.path.expanduser(loc_csv_path)
    output_path = os.path.expanduser(output_path)
    os.makedirs(output_path, exist_ok=True)

    # Load loc metadata
    loc_df = pd.read_csv(loc_csv_path)
    loc_df['longitude_wgs'] = loc_df['longitude_wgs'].astype(float)
    loc_df['latitude_wgs'] = loc_df['latitude_wgs'].astype(float)
    if loc_id_list is not None:
        loc_df = loc_df[loc_df['loc_id'].isin(loc_id_list)]

    # For each year and NetCDF file
    for year, file_path in nc_file_paths.items():
        file_path = os.path.expanduser(file_path)
        print(f"Processing ERA5 data for {year}...")

        try:
            ds = Dataset(file_path)
        except Exception as e:
            print(f"[Error] Cannot open {file_path}: {e}")
            continue

        # Determine variables to extract
        if variables == 'all':
            # Exclude coordinate and time variables
            variables_to_extract = [v for v in ds.variables if v not in ['lat', 'lon', 'time']]
            print(f"[Auto] Extracting variables: {variables_to_extract}")
        else:
            variables_to_extract = variables

        time_var = ds.variables['valid_time']
        dates = pd.to_datetime(nc4num2date(time_var[:], time_var.units))
        lon = np.array(ds['longitude'][:])
        lat = np.array(ds['latitude'][:])

        for loc_id, row in loc_df.set_index('loc_id').iterrows():
            loc_file = os.path.join(output_path, f"loc_id_{loc_id}_{year}.csv")
            if not os.path.exists(loc_file):
                print(f"[Skip] No pollutant file for loc_id {loc_id}, skipping...")
                continue

            loc_data = pd.read_csv(loc_file)
            lon0 = row['longitude_wgs']
            lat0 = row['latitude_wgs']
            lon_idx = np.abs(lon - lon0).argmin()
            lat_idx = np.abs(lat - lat0).argmin()

            for var in variables_to_extract:
                print(f"  Extracting {var} for loc {loc_id} ({year})")
                try:
                    var_data = np.array(ds[var][:, lat_idx, lon_idx], dtype=np.float32)
                except Exception as e:
                    print(f"    [Warning] Variable {var} not found or error reading: {e}")
                    continue

                
                if len(loc_data) != len(var_data):
                    print(f"    [Mismatch] {len(var_data)} values for {var}, but {len(loc_data)} rows. Skipping.")
                    continue

                
                loc_data[var] = var_data

            # Save back updated CSV
            loc_data.to_csv(loc_file, index=False)

        ds.close()

def nc4num2date(times, units):
    from netCDF4 import num2date
    return num2date(times, units).astype(str)
def extract_loc_all_data(
    loc_id,
    pollutants='all',
    pollutant_years=['2023', '2024'],
    meteorology_years=['2023', '2024'],
    nc_pollutant_root="/rds/projects/l/liub-eulez/City",
    nc_meteorology_paths=None,
    loc_csv_path="/rds/projects/l/liub-eulez/City/loc_info.csv",
    output_path="/rds/projects/l/liub-eulez/City/Results"
):
    """
    Extract all data (pollutants + meteorology) for a given loc_id.

    Parameters:
        loc_id (int): Loc point ID to extract
        pollutants (list or 'all'): List of pollutant names or 'all'
        pollutant_years (list): Years of pollutant data
        meteorology_years (list): Years of meteorology data
        nc_pollutant_root (str): Folder containing pollutant-year folders
        nc_meteorology_paths (dict): Mapping from year to nc path for meteorology
        loc_csv_path (str): Path to loc info CSV
        output_path (str): Output folder for loc_id CSV
    """
    import pandas as pd
    import os

   
    extract_pollutants(
        nc_root_path=nc_pollutant_root,
        loc_csv_path=loc_csv_path,
        pollutants=pollutants,
        years=pollutant_years,
        output_path=output_path,
        loc_id_list=[loc_id]
    )

   
    if nc_meteorology_paths is None:
        nc_meteorology_paths = {
            y: os.path.join(nc_pollutant_root, "ERA5City_Merged", f"ERA5{y}_daily.nc")
            for y in meteorology_years
        }

    extract_meteorology(
        nc_file_paths=nc_meteorology_paths,
        loc_csv_path=loc_csv_path,
        variables='all',
        output_path=output_path,
        loc_id_list=[loc_id]
    )

  
    print(f"[Merge] Combining years for loc_id {loc_id}...")
    dfs = []
    for year in sorted(set(pollutant_years + meteorology_years)):
        file_path = os.path.join(os.path.expanduser(output_path), f"loc_id_{loc_id}_{year}.csv")
        if os.path.exists(file_path):
            df = pd.read_csv(file_path, parse_dates=['date'])
            dfs.append(df)
        else:
            print(f"[Skip] File not found: {file_path}")

    if dfs:
        merged = pd.concat(dfs, ignore_index=True).sort_values("date")
        merged_file = os.path.join(os.path.expanduser(output_path), f"loc_id_{loc_id}.csv")
        merged.to_csv(merged_file, index=False)
        print(f"[All Done] Merged and saved: {merged_file}")
    else:
        print(f"[Error] No yearly files found to merge for loc_id {loc_id}")
#%%
#extract_loc_all_data(
#    loc_id=6,
#    pollutant_years=['2023'],
#    meteorology_years=['2023'],
#    nc_pollutant_root="/path/to/your/Downloads/City",
#    loc_csv_path="/path/to/your/Downloads/City/loc_info.csv",
#    output_path="/path/to/your/Downloads/City/results"
#)