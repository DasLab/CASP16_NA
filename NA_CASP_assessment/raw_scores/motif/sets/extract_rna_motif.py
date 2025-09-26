import pandas as pd
from glob import glob
import os
from tqdm import tqdm

def parse_rna_motif(input_directory, name):
    # Initialize a list to hold our data
    data,errors,no_motifs,incomplete_runs = [],[],[],[]

    # Iterate over all stdout files
    for file_path in tqdm(glob(input_directory)):
        file_name = os.path.basename(file_path)
        stderr_file_path = file_path.replace("stdout", "stderr")  
        is_error = False
        if os.path.getsize(stderr_file_path) > 0:
            is_error = True
            # The stderr file is not empty, log the error file
            # now only do if not motifs dounf ierrors.append({'file_name': file_name, 'error_file': stderr_file_path})
            #continue  # Skip to the next stdout file
        with open(file_path, 'r') as file:
            motifs_found = False
            incomplete_run = True
            for line in file:
                if "rna_motif: Created: color_motifs.pml" in line:
                    incomplete_run = False
                # Check if line starts with spaces
                if line.startswith(" "):
                    motifs_found = True
                    parts = line.strip().split()
                
                    # The first part is the motif type (remove the last character, which is ":")
                    motif_type = parts[0][:-1]  # Remove the colon
                

                    # Process residues, with the edge cases that are qeirdly fromatting, guessing segements
                    residues = []
                    for i in range(1, len(parts)):
                        if parts[i].startswith(':') and i > 0:  # Check if it is a special base
                            # Combine with the previous residue
                            combined = f"{residues[-1]}{parts[i]}"
                            residues[-1] = combined
                        else:
                            residues.append(parts[i])
                    # The remaining parts are joined with "_"
                    residues = '_'.join(residues).replace("[RAD:deoxy_O2prime]", "").replace("[RCY:deoxy_O2prime]","")  # Join the residues parts

                    # Append the data to our list
                    data.append({'file_name': file_name, 'motif_type': motif_type, 'residues': residues})
            if not motifs_found:
                no_motifs.append({'file_name': file_name})
                if is_error:
                    errors.append({'file_name': file_name, 'error_file': stderr_file_path})
            if incomplete_run:
                incomplete_runs.append({'file_name': file_name})

    # Create pandas DataFrames from the data lists
    df_motifs = pd.DataFrame(data)
    df_no_motifs = pd.DataFrame(no_motifs)
    df_errors = pd.DataFrame(errors)
    df_incomplete_runs = pd.DataFrame(incomplete_runs)

    # Save the DataFrames to CSV files
    df_motifs.to_csv(f'{name}_motif_data.csv', index=False)
    df_no_motifs.to_csv(f'{name}_no_motifs.csv', index=False)
    df_errors.to_csv(f'{name}_errors.csv', index=False)
    df_incomplete_runs.to_csv(f'{name}_incomplete_runs.csv', index=False)

parse_rna_motif("rna_motif_outputs/RNA_Data_rna_motif/*stdout","predicted")
parse_rna_motif("rna_motif_outputs/Targets_RNA_rna_motif/*stdout","native")

