########################################

# The code below downloads and processes ITS_LIVE data for specified glacier coordinates.
# It reads coordinates from CSV files, finds the corresponding granule URLs, fetches timeseries velocity data for a given latitude and longitude,
# Finally, it exports the results to new CSV files in a structured output directory.


# The code uses the ITS_LIVE (https://its-live.jpl.nasa.gov) Python library to access the data cubes and perform the necessary operations.
# This file uses functions from the ITS_LIVE project (https://github.com/nasa/itslive), licensed under GNU GPLv3.

########################################

import os
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import sys
import re
import xarray as xr

sys.path.append(os.path.join(os.path.dirname(__file__), 'itslive-py'))

from itslive.velocity_cubes._cubes import get_time_series, export_csv

# Define directories and catalog URL
CATALOG_URL = "https://its-live-data.s3.amazonaws.com/datacubes/catalog_v02.json"
INPUT_DIR = "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/New Points v3/Points/Sampling Points/as CSV"
OUTPUT_DIR = "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Jess ITS_LIVE v3"

# Load ITS_LIVE catalog
catalog = gpd.read_file(CATALOG_URL)

# Function to find the granule URL for given latitude and longitude
def find_granule_url(lat, lon):
    point = gpd.GeoSeries([Point(lon, lat)], crs="EPSG:4326")
    gdf = gpd.GeoDataFrame(geometry=point)
    match = catalog.sjoin(gdf, how="inner", predicate="intersects")
    return match.iloc[0]['zarr_url'] if not match.empty else None

# Function to extract glacier information from path
def extract_glacier_info(path):
    filename = os.path.basename(path)
    match = re.match(r"gl_(\d+)_(\d+)_(\d+)\.csv", filename)
    if match:
        glacier, flowline, elevation = match.groups()
        return f"Glacier {glacier}, Flowline {flowline}, Elevation {elevation}"
    return "Unknown Glacier"

# Function to fetch timeseries velocity data and export to CSV
def process_csv(input_csv_path, output_csv_path):
    glacier_info = extract_glacier_info(input_csv_path)
    df = pd.read_csv(input_csv_path)
    all_points = []

    for _, row in df.iterrows():
        lat, lon = row['y_4326'], row['x_4326']
        url = find_granule_url(lat, lon)
        if not url:
            print(f"!!! No granule found for ({lat}, {lon}) in {glacier_info}")
            continue

        try:
            result_series = get_time_series([(lon, lat)])
            if not result_series:
                print(f"!!! Empty time series for ({lat}, {lon} in {glacier_info})")
                continue
            all_points.append((lon, lat))
        except Exception as e:
            print(f"!!! Error fetching data for ({lat}, {lon}) in {glacier_info}: {e}")

    if all_points:
        export_csv(all_points, outdir=output_csv_path)
    else:
        print(f"!!! No data exported for {glacier_info}")

# Main processing loop
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for glacier_folder in os.listdir(INPUT_DIR):
        g_path = os.path.join(INPUT_DIR, glacier_folder)
        if not os.path.isdir(g_path):
            continue

        for file in os.listdir(g_path):
            if not file.endswith(".csv"):
                continue

            input_csv = os.path.join(g_path, file)
            name_out = file.replace(".csv", "_itslive_v3")
            output_dir_path = os.path.join(OUTPUT_DIR, glacier_folder, name_out)
            os.makedirs(output_dir_path, exist_ok=True)
            glacier_info = extract_glacier_info(input_csv)

            print(f"Processing: {glacier_info}")
            process_csv(input_csv, output_dir_path)

if __name__ == "__main__":
    main()