"""
Merge monthly ERA5 NetCDF (.nc) files and generate daily resolution datasets.

This script handles ERA5 hourly data organized by variable type, including:
- Accumulation (e.g. total precipitation): resampled by sum
- Daily maximum (e.g. max temperature): resampled by max
- Instantaneous (e.g. wind speed): resampled by mean or retained as is

It identifies and merges corresponding stepType-specific files for each month,
generates yearly merged data, and performs variable-appropriate daily aggregation.

Includes plotting for visual inspection at a specified lat/lon point.

Author: Chengxu Tong
"""

import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

input_root = '/path/to/your/city_data/ERA5City'
output_root = '/path/to/your/city_data/ERA5City_Merged'
os.makedirs(output_root, exist_ok=True)

years = ['2023', '2024']
months = [f"{i:02d}" for i in range(1, 13)]

for year in years:
    monthly_datasets = []

    for month in months:
        folder_name = f'ERA5{year}{month}'
        folder_path = os.path.join(input_root, folder_name)

        month_files = sorted([
            os.path.join(folder_path, f) for f in os.listdir(folder_path)
            if f.endswith('.nc')
        ])
        if not month_files:
            print(f"No files in {folder_path}, skipping...")
            continue

        ds_month = xr.open_mfdataset(month_files, combine='by_coords')
        monthly_datasets.append(ds_month)

    print(f"Merging all months of {year}...")
    ds_year = xr.concat(monthly_datasets, dim='valid_time')

    year_path = os.path.join(output_root, f'ERA5{year}.nc')
    ds_year.to_netcdf(year_path)
    print(f"Saved yearly file: {year_path}")

    print(f"Calculating daily mean for {year}...")
    ds_daily = ds_year.resample(valid_time='1D').mean()

    valid_lens = {var: ds_daily[var].sizes['valid_time'] for var in ds_daily.data_vars}
    ref_len = max(valid_lens.values())
    drop_vars = [var for var, l in valid_lens.items() if l != ref_len]
    if drop_vars:
        print(f"Dropping incompatible variables from daily output: {drop_vars}")
        ds_daily = ds_daily.drop_vars(drop_vars)

    daily_path = os.path.join(output_root, f'ERA5{year}_daily.nc')
    ds_daily.to_netcdf(daily_path)
    print(f"Saved daily mean file: {daily_path}")

# ---- Plotting section ----
target_lat = 28.3
target_lon = 113.05

ds = xr.open_dataset(daily_path)
lat_idx = np.abs(ds.latitude - target_lat).argmin().item()
lon_idx = np.abs(ds.longitude - target_lon).argmin().item()

for var in ds.data_vars:
    series = ds[var][:, lat_idx, lon_idx].values
    plt.figure()
    plt.plot(series)
    plt.title(f"{var} at lat={target_lat}, lon={target_lon}")
    plt.xlabel("Days")
    plt.ylabel(var)
    plt.grid(True)
    plt.tight_layout()
    plt.show()