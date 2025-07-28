########################################

# The code below downloads and processes ITS_LIVE data for specified glacier coordinates.
# It reads coordinates from CSV files, finds the corresponding granule URLs, fetches time series data,
# and exports the results to new CSV files named after the glacier folder with "_itslive_v2.csv" appended to the name.


# The code uses the ITS_LIVE (https://its-live.jpl.nasa.gov) Python library to access the data cubes and perform the necessary operations.
# This file uses functions from the ITS_LIVE project (https://github.com/nasa/itslive), licensed under GNU GPLv3.

########################################

import os
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import sys

# Make sure we can import from local itslive-py
sys.path.append(os.path.join(os.path.dirname(__file__), 'itslive-py'))

from itslive.velocity_cubes._cubes import get_time_series, export_csv

# Constants
CATALOG_URL = "https://its-live-data.s3.amazonaws.com/datacubes/catalog_v02.json"
INPUT_DIR = "Coordinates"
OUTPUT_DIR = "Output"

# Load ITS_LIVE catalog
catalog = gpd.read_file(CATALOG_URL)

def find_granule_url(lat, lon):
    point = gpd.GeoSeries([Point(lon, lat)], crs="EPSG:4326")
    gdf = gpd.GeoDataFrame(geometry=point)
    match = catalog.sjoin(gdf, how="inner", predicate="intersects")
    return match.iloc[0]['zarr_url'] if not match.empty else None

def process_csv(input_csv_path, output_csv_path):
    df = pd.read_csv(input_csv_path, header=None, names=["latitude", "longitude"])
    all_data = []

    for _, row in df.iterrows():
        lat, lon = row['latitude'], row['longitude']
        url = find_granule_url(lat, lon)
        if not url:
            print(f"⚠️ No granule found for ({lat}, {lon})")
            continue

        try:
            result = get_time_series(lat, lon)
            if result is None or result.dims.get("time", 0) == 0:
                print(f"⚠️ Empty time series for ({lat}, {lon})")
                continue

            result.attrs["lat"] = lat
            result.attrs["lon"] = lon
            all_data.append(result)

        except Exception as e:
            print(f"Error fetching data for ({lat}, {lon}): {e}")

    if all_data:
        export_csv(all_data, output_csv_path)
    else:
        print(f"No data exported for {input_csv_path}")

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for glacier_folder in os.listdir(INPUT_DIR):
        g_path = os.path.join(INPUT_DIR, glacier_folder)
        if not os.path.isdir(g_path):
            continue

        for file in os.listdir(g_path):
            if not file.endswith("_coord.csv"):
                continue

            input_csv = os.path.join(g_path, file)
            name_out = file.replace("_coord.csv", "_itslive_v2.csv")
            glacier_output_path = os.path.join(OUTPUT_DIR, glacier_folder)
            os.makedirs(glacier_output_path, exist_ok=True)
            output_csv = os.path.join(glacier_output_path, name_out)

            print(f"Processing: {input_csv}")
            process_csv(input_csv, output_csv)
            print(f"Saved: {output_csv}")

if __name__ == "__main__":
    main()