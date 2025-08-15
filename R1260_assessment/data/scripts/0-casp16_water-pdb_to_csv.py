from glob import glob
from os import rename
import argparse
import numpy as np
from tqdm import tqdm
import pandas as pd
from biopandas.pdb import PandasPdb

'''
python scripts/0-casp16_water-pdb_to_csv_single_file.py \
    --input_folder R1260 \
    --output_folder R1260_csvs_simple \
    --min_dist 5 \
    --csv_output summary_atom_counts_simple.csv \
    --group 156
'''

parser = argparse.ArgumentParser(description="Process PDB files and generate summary and CSV outputs.")
parser.add_argument('--input_folder', required=True, help="Folder containing the folders of PDB files.")
parser.add_argument('--output_folder', required=True, help="Path to save the summary output CSV file.")
parser.add_argument('--csv_output', required=True, help="Foldser to save the individual CSV files.")
parser.add_argument('--min_dist', type=float, default=5, help="Cutoff to ignore solvent (angstrom).")
parser.add_argument('--group', default=None, type=str, help="If specified only process files from that group.")
args = parser.parse_args()

desired_columns = ['record_name', 'atom_name',  'residue_name', 'residue_number', 'x_coord', 'y_coord', 'z_coord', 'element_symbol']
residue_dict = {'A':'A', 'C':'C', 'G':'G', 'U':'U',
                'RA':'A', 'RC':'C', 'RG':'G', 'RU':'U',
                'G5':'G', 'G3':'G', 'RG5':'G', 'RG3':'G',
                'HOH':'HOH', 'WAT':'HOH', 'SOL':'HOH',
                'MG':'MG', 'mMg':'MG','nMg':'MG',
                'NA':'NA', 'Na+':'NA',
                'CL':'CL', 'Cl-':'CL','UNL':'CL', 
                'K':'K'}
RNA_residues = ['A','C','G','U']

def dist_minRef_to_all(ref_coords, coords):
    dists = np.linalg.norm(ref_coords[:, np.newaxis, :] - coords[np.newaxis, :, :], axis=-1)
    min_dists = dists.min(axis=0)
    return min_dists

data = {'file_location':[],'num_H':[],'atoms_A':[],'atoms_C':[],'atoms_G':[],'atoms_U':[]}
current_residue_list = ['A','C','G','U']

# this assume files are already untarred
if args.group:
    pdb_folders = glob(f"{args.input_folder}/R1260TS{args.group}")
else:
    pdb_folders = glob(f"{args.input_folder}/R1260TS???")
for pdb_folder in pdb_folders:
    
    pdb_files = glob(f"{pdb_folder}/*pdb")

    # add the pdb suffix when necessary
    if len(pdb_files) == 0:
        pdb_files = glob(f"{pdb_folder}/R1260TS???_*")
        for pdb_file in pdb_files:
            rename(pdb_file, f"{pdb_file}.pdb")
        pdb_files = glob(f"{pdb_folder}/*pdb")
    
    print(f"Working on {len(pdb_files)} structures from {pdb_folder}.")
    pdb_files = sorted(pdb_files, key=lambda file_name: int(file_name.rsplit('_', 1)[1].split('.')[0]))

    rna_dfs = []
    sol_dfs = []
    for pdb_file in tqdm(pdb_files):
        # open their pdb file
        data['file_location'].append(pdb_file)
        ppdb = PandasPdb().read_pdb(pdb_file)
        
        # just assume their ATOM HETATM differentiation means nothing
        df = pd.concat([ppdb.df['ATOM'][desired_columns],ppdb.df['HETATM'][desired_columns]])

        # Check for unknown residue names
        unknown_residues = df.loc[~df['residue_name'].isin(residue_dict.keys()), 'residue_name'].unique()
        if len(unknown_residues) > 0:
            raise ValueError(f"Unknown residue names found: {unknown_residues}")

        # Rename to standard names
        df['residue_name'] = df.residue_name.replace(residue_dict)

        # remove H
        data['num_H'].append((df.element_symbol=='H').sum())
        df = df[(df.element_symbol != 'H') & (df.atom_name.str[0]!='H')]

        # check what other residues there are
        residues = list(df.residue_name.unique())
        non_RNA_residues = [x for x in residues if x not in RNA_residues]

        # remove up to 5A from any RNA molecule
        rna_df = df[df.residue_name.isin(RNA_residues)].copy()
        sol_df = df[df.residue_name.isin(non_RNA_residues)].copy()
        rna_coords = rna_df[['x_coord','y_coord','z_coord']].to_numpy().astype(float)
        sol_coords = sol_df[['x_coord','y_coord','z_coord']].to_numpy().astype(float)
        sol_dist = dist_minRef_to_all(rna_coords,sol_coords)
        sol_close_df = sol_df[sol_dist <= args.min_dist]

        # save csvs
        # rna_df.to_csv(f'{args.output_folder}/rna_{pdb_file.rsplit("/",1)[1]}', index=False)
        # sol_df.to_csv(f'{args.output_folder}/sol_{pdb_file.rsplit("/",1)[1]}', index=False)
        frame = pdb_file.rsplit("/",1)[1].split('_')[-1].rsplit('.')[0]
        rna_df['frame'] =  frame
        sol_df['frame'] = frame
        rna_dfs.append(rna_df)
        sol_dfs.append(sol_df)

        number_of_atoms = df.groupby('residue_name').size().to_dict()
        for residue in set(residues+current_residue_list):
            if f"atoms_{residue}" in data:
                if residue in number_of_atoms:
                    data[f"atoms_{residue}"].append(number_of_atoms[residue])
                else:
                    data[f"atoms_{residue}"].append(0)
            else:
                current_residue_list.append(residue)
                data[f"atoms_{residue}"] = [0]*(len(data['num_H'])-1)+[number_of_atoms[residue]]
    group_name = pdb_file.rsplit("/",1)[1].split('_')[0]
    pd.concat(rna_dfs).to_csv(f'{args.output_folder}/rna_{group_name}.csv', index=False)
    pd.concat(sol_dfs).to_csv(f'{args.output_folder}/sol_{group_name}.csv', index=False)

df = pd.DataFrame(data)
df['group'] = df.file_location.apply(lambda x: x.rsplit('_',1)[0][-3:])
df['model'] = df.file_location.apply(lambda x: x.rsplit('_',1)[1][:-4])
df = df.astype({'model':int})
df = df.sort_values(by=['group','model'])
df.to_csv(args.csv_output, index=False)
