from pathlib import Path

base_dir = Path("/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Jess ITS_LIVE v2/Inputs/Observations/v2")

for glacier_folder in base_dir.glob("Gl*"):
    if not glacier_folder.is_dir():
        continue

    for band_folder in glacier_folder.glob("Gl*_v2"):
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

        print(f"✅ Combined to: {combined_path}")