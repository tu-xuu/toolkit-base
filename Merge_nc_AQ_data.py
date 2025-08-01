#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility script for merging daily pollutant .nc files into monthly and yearly datasets.
Suitable for processing Jing Wei's High Air Pollutants daily NetCDF data (e.g., NO2, PM2.5).
Author: Chengxu Tong
"""
import os
import xarray as xr
import pandas as pd
import re
from tqdm import tqdm

def extract_date_from_filename(filename):
    """
    Extracts date from filename in format "_D1K_YYYYMMDD_".
    Returns a pandas.Timestamp object if matched, otherwise None.
    """
    match = re.search(r'_D1K_(\d{8})_', filename)
    if match:
        return pd.to_datetime(match.group(1), format="%Y%m%d")
    return None

def merge_one_pollutant_year_safe(pollutant, year, base_dir="/path/to/your/city_data"):
    """
    Merge .nc files for a given pollutant and year into monthly NetCDF files, then merge monthly files.
    - Reads files from a pollutant-year-specific folder.
    - Extracts date from filenames and groups by month.
    - Converts variables to float32 for memory efficiency.
    - Saves monthly files to a temporary folder, then merges them.
    Parameters:
        pollutant (str): Pollutant name (e.g., "NO2")
        year (int): Year to process (e.g., 2022)
        base_dir (str): Base directory containing pollutant data
    """
    input_dir = os.path.join(base_dir, f"{pollutant}_{year}")
    temp_dir = os.path.join(base_dir, f"{pollutant}_temp_{year}")
    output_dir = os.path.join(base_dir, f"{pollutant}merge")
    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    file_list = sorted([
        f for f in os.listdir(input_dir)
        if f.endswith(".nc") and f"{year}" in f
    ])

    month_groups = {}
    for fname in file_list:
        date = extract_date_from_filename(fname)
        if date is not None:
            month_key = date.strftime("%Y-%m")
            month_groups.setdefault(month_key, []).append(fname)

    
    monthly_files = []
    for month, files in tqdm(month_groups.items(), desc=f"{pollutant}_{year}"):
        datasets = []
        for fname in sorted(files):
            fpath = os.path.join(input_dir, fname)
            ds = xr.open_dataset(fpath)
            date = extract_date_from_filename(fname)
            if date is not None:
                ds = ds.expand_dims("time")
                ds["time"] = [date]
                for var in ds.data_vars:
                datasets.append(ds)
        if datasets:
            monthly_ds = xr.concat(datasets, dim="time")
            monthly_path = os.path.join(temp_dir, f"{month}.nc")
            monthly_ds.to_netcdf(monthly_path, format="NETCDF3_64BIT")
            monthly_files.append(monthly_path)

   
    if monthly_files:
        all_months = [xr.open_dataset(f) for f in monthly_files]
        merged = xr.concat(all_months, dim="time")
        output_path = os.path.join(output_dir, f"{pollutant}_{year}_merged.nc")
        merged.to_netcdf(output_path, format="NETCDF3_64BIT")
        print(f"✅ Final saved: {output_path}")


merge_one_pollutant_year_safe("CO", "2023")