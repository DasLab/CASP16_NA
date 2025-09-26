import pandas as pd
import glob

def filter_canonical_bps(input_pattern="base_pair_tables/*_combined_dssr_data.csv", output_suffix="_canonical"):

    for csv_file in glob.glob(input_pattern):
        df = pd.read_csv(csv_file)
        df_filtered = df[df['name'].isin(['WC', 'Wobble'])]  # efficient filtering

        output_filename = csv_file.replace(".csv", f"{output_suffix}.csv")
        df_filtered.to_csv(output_filename, index=False)
        
        df_non = df[~df.name.isin(['WC','Wobble'])]
        df_non.to_csv(csv_file.replace(".csv",f"not{output_suffix}.csv"),index=False)
if __name__ == "__main__":
    filter_canonical_bps()
