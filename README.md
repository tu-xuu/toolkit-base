# toolkit-base

**Personal code snippets for data analysis and modeling**

---

## Extracting Pollutant and Meteorological Data from NetCDF file

This section contains scripts to extract and process daily air pollutant and meteorological data 
from NetCDF files at specified latitude/longitude locations.

### Features

- Extract pollutants (e.g., NO2, PM2.5) and meteorology by location
- Merge yearly data for each location into CSV format
- Parallel processing for large-scale extraction
- Optional visualization of trends at selected points

### File Overview

| File | Description |
|------|-------------|
| [`Extract_Pollutant_Met_Location.py`](./Extract_Pollutant_Met_Location.py) | Main extractor function script |
| [`ERA5_Extract_DailyAvg_Plot.py`](./ERA5_Extract_DailyAvg_Plot.py)         | ERA5 daily merge and plot |
| [`run_parallel_extract.py`](./run_parallel_extract.py)           | Multiprocessing runner |
| [`Merge_nc_AQ_data.py`](./Merge_nc_AQ_data.py)               | Optional: Merge daliy AQ data to monthly or yearly|

---

## Other Tools Coming Soon

This repository will continue to grow with additional utilities for data cleaning, 
analysis, modeling, and visualization.
