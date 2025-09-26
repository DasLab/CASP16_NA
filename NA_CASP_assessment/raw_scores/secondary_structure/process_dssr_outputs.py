import os
import subprocess
import glob
import pandas as pd

def process_dssr_outputs(dssr_output_dir="dssr_outputs",dssr_out_file="output.txt"):
    """
    Processes DSSR output files in subfolders, handles missing files, and combines CSV results.
    """

    csv_files = []
    no_output_folders = []

    for root, _, files in os.walk(dssr_output_dir):
        output_txt_path = os.path.join(root, dssr_out_file) 
        if os.path.exists(output_txt_path):
            # Check if output.txt is empty
            if os.stat(output_txt_path).st_size == 0:
                no_output_folders.append(root)
                continue # Skip to the next folder

            try:
                # run the bash script to get base pairs
                command = ["./parse_dssr_bp_list.sh", output_txt_path, os.path.join(root, "output.csv")]
                subprocess.run(command, check=True, capture_output=True, text=True)
                csv_files.append(os.path.join(root, "output.csv"))
            except subprocess.CalledProcessError as e:
                print(f"Error processing {root}: {e}")  # Handle errors during shell execution
            except Exception as e:
                print(f"An unexpected error occurred processing {root}: {e}")
        else:
            no_output_folders.append(root)

    # Write list of folders without or with empty output.txt to a file.
    with open("no_dssr_output.txt", "w") as f:
        for folder in no_output_folders:
            f.write(folder + "\n")

    # get all targets to save a csv per target
    all_folders = [f.path for f in os.scandir(dssr_output_dir) if f.is_dir()]
    r_folders = [os.path.basename(f) for f in all_folders if os.path.basename(f).startswith("R")]
    targets = set()
    for folder in r_folders:
        target = folder.split("TS")[0].split("LG")[0]
        targets.add(target)
    print("Targets:", targets)

    # combine csvs
    for target in targets:
       combined_data = []
       for csv_file in csv_files:
            if target in csv_file:
                try:
                    df = pd.read_csv(csv_file)
                    folder_name = os.path.basename(os.path.dirname(csv_file))
                    df['folder'] = folder_name  # Add folder name as a new column
                    combined_data.append(df)
                except pd.errors.EmptyDataError:
                    print(f"Warning: CSV file {csv_file} is empty and will be skipped.")
                except pd.errors.ParserError:
                    print(f"Warning: Could not parse CSV file {csv_file} and will be skipped.")
                except Exception as e:
                    print(f"An unexpected error occurred reading {csv_file}: {e}")
       combined_df = pd.concat(combined_data, ignore_index=True)
       combined_df.to_csv(f"base_pair_tables/{target}_combined_dssr_data.csv", index=False)

if __name__ == "__main__":
    process_dssr_outputs()
