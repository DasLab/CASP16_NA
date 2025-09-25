import os
import pandas as pd
import re


def extract_interaction_info(filename):
    """Extracts interaction data, handling motifs files carefully."""
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"File not found: {filename}")
        return None

    interaction_type = os.path.basename(filename).split(".")[2]

    interactions = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) < 2:  # Skip lines with insufficient data
            continue

        residues = []
        chains = []
        if interaction_type == 'motifs':
            
            for part in parts[1:]:
                chains.append(part.split(':')[0])
                if "-" in part:  # Handle ranges like 481-483
                    try:
                        start, end = map(int, part.split(":")[1].split("-"))
                        residues.extend(range(start, end + 1))
                    except (ValueError, IndexError):
                        print(f"Warning: Skipping invalid residue range: {part} in file {filename}")
                        continue
                elif ":" in part:  # Handle single residues like 0:485
                    try:
                        residues.append(int(part.split(":")[1]))
                    except (ValueError, IndexError):
                        print(f"Warning: Skipping invalid residue: {part} in file {filename}")
                        continue
                else:
                    print(f"Warning: Skipping unexpected part: {part} in file {filename}")
            if "0" in chains and "1" in chains and all(29 <= res <= 371 for res in residues): # any
                interactions.append(" ".join(parts))  # Include interaction type
        else:
            residues.append(int(parts[0].split(":")[1]))
            residues.append(int(parts[1].split(":")[1]))
            chains.append(int(parts[0].split(":")[0]))
            chains.append(int(parts[1].split(":")[0]))
            if 0 in chains and 1 in chains and all(29 <= res <= 371 for res in residues):
                interactions.append(" ".join(parts))  # Include interaction type

    return interaction_type, interactions


def process_prediction_files(directory):
    """Processes all relevant files in the directory."""
    data = []
    for filename in os.listdir(directory):
        if filename.startswith('R1285TS') and filename.endswith(('.base_pairs.txt','.other_contacts.txt','.stacks.txt','.motifs.txt')): #Added more file extensions
            filepath = os.path.join(directory, filename)
            prediction_name = filename.split(".")[0]  #Removes the file extensions
            interaction_info = extract_interaction_info(filepath)
            if interaction_info:
                interaction_type, interactions = interaction_info
                for interaction in interactions:
                    data.append([prediction_name, interaction_type, interaction])

    return pd.DataFrame(data, columns=['prediction', 'interaction_type', 'interaction_details'])

directory_path = "R1285o"  
df = process_prediction_files(directory_path)
print(df)
df.to_csv("all_intermolecular_interactions.csv", index=False)
