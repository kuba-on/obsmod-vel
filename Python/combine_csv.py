########################################

# The code below combines multiple CSV files from each glacier's band folder into a single CSV file.
# It reads all CSV files in each band folder, skips headers after the first file, and writes the combined data to a new CSV file named after the band folder with "_comb.csv" appended to the name.
# The combined CSV files are saved in the respective glacier folder.

########################################

from pathlib import Path

base_dir = Path("/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Inputs/Observations")

for glacier_folder in base_dir.glob("Glacier*"):
    if not glacier_folder.is_dir():
        continue

    for band_folder in glacier_folder.glob("gl*_v3"):
        if not band_folder.is_dir():
            continue

        combined_name = band_folder.name + "_comb.csv"
        combined_path = glacier_folder / combined_name

        with combined_path.open("w") as outfile:
            first = True
            for csv_file in sorted(band_folder.glob("*.csv")):
                with csv_file.open("r") as infile:
                    lines = infile.readlines()
                    if not lines:
                        continue
                    if first:
                        outfile.writelines(lines)
                        first = False
                    else:
                        outfile.writelines(lines[1:])  # skip header

        print(f" Combined to: {combined_path}")