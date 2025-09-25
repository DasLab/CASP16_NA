import os
import re
import csv

def process_interaction_line(line):
    parts = line.split()
    chain1, residue1 = parts[0].split(":")
    chain2, residue2 = parts[1].split(":")
    interaction_codes = parts[2:]  #Get interaction codes
    return chain1, residue1, chain2, residue2, interaction_codes


def extract_interaction_data(filename, interaction_type): #interaction_type added as argument
    interactions = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            #chain1, residue1, chain2, residue2, interaction_codes = process_interaction_line(line)
            #interactions.append(f"{chain1}:{residue1} {chain2}:{residue2} {' '.join(interaction_codes)}") #Reconstruct interaction string
            interactions.append(line)
    return interactions


def rename_chains(interaction_details):
    return interaction_details.replace("A:", "0:").replace("B:", "1:")


native_files = ["native/B1_native.txt", "native/B2_native.txt", "native/B3_native.txt","native/kink_native.txt"] # 
r1285_dir = "R1285o"

results = []

for native_file in native_files:
    motif = native_file.split("_")[0]
    with open(native_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip()
            if not line:
                continue

            parts = line.split(":")
            if len(parts) < 2:
                continue
            
            # OLE.pdb.base_pairs.txt:A:315 B:318 W W C

            native_interaction_type = parts[0].split(".")[2]
            native_interaction_details = ':'.join(parts[1:])
            for filename in os.listdir(r1285_dir):
                interaction_type_from_filename = filename.split(".")[-2]

                if interaction_type_from_filename == native_interaction_type: # Match Interaction Type
                    base_filename = filename[:-len(f".{native_interaction_type}.txt")]

                    filepath = os.path.join(r1285_dir, filename)
                    #try:
                    r1285_interactions = extract_interaction_data(filepath,native_interaction_type)
                    renamed_native_interaction = rename_chains(native_interaction_details)

                    present = renamed_native_interaction in r1285_interactions
                    results.append({
                        "prediction": base_filename,
                        "motif": motif,
                        "interaction_type": native_interaction_type,
                        "residue": renamed_native_interaction,#native_interaction_details,
                        "present": present
                    })

# Save to CSV
with open('interaction_results.csv', 'w', newline='') as csvfile:
    fieldnames = ['prediction', 'motif', 'interaction_type', 'residue', 'present']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)

print("Results saved to interaction_results.csv")