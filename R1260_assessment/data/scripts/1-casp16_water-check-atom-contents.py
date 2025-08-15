from biopandas.pdb import PandasPdb
import pandas as pd
import argparse
from tqdm import tqdm
from glob import glob
import os

'''
python scripts/1-casp16_water-check-atom-contents.py \
    --native_pdb 9CBU_cleaned.pdb \
    --input_folder R1260_serial \
    --native_output scripts/9CBU-rna.csv
'''

# for the simplified RNA pdb, check all atoms present
parser = argparse.ArgumentParser(description="Check CSV outputs.")
parser.add_argument('--native_pdb', required=True, help="The native pdb containing rna to align to aka reference.")
parser.add_argument('--native_output', required=True, help="The reference csv to output.")
parser.add_argument('--input_folder', required=True, help="Folder of the predicted CSV files.")
args = parser.parse_args()

# expected atoms
backbone_atoms = ['P','OP1','OP2',"O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'"]
expected_atoms = {'G':backbone_atoms+['N9','C8','N7','C5','C6','O6','N1','C2','N2','N3','C4'],
                  'U':backbone_atoms+['N1','C2','O2','N3','C4','O4','C5','C6'],
                  'C':backbone_atoms+['N1','C2','O2','N3','C4','N4','C5','C6'],
                  'A':backbone_atoms+['N9','C8','N7','C5','C6','N6','N1','C2','N3','C4']}
# G3 normal, if first residue, no OP1, OP2, P

# Reference
# check cleaned pdb and order
# convert to csv
desired_columns = ['record_name', 'atom_name',  'residue_name', 'residue_number', 'x_coord', 'y_coord', 'z_coord', 'element_symbol']
RNA_residues = ['A','C','G','U']

ppdb = PandasPdb().read_pdb(args.native_pdb)
df = pd.concat([ppdb.df['ATOM'][desired_columns],ppdb.df['HETATM'][desired_columns]])

# remove H
df = df[df.element_symbol != 'H']

# order
df = df[~((df.residue_number==22)&(df.atom_name.isin(["O5'",'OP1','OP2','P'])))]
df = df.sort_values(by=['residue_number','atom_name'])

# save csvs
df[df.residue_name.isin(RNA_residues)].to_csv(args.native_output, index=False)

df = pd.read_csv(args.native_output)
first_res = df.residue_number.min()

# check content and order
for resn,resn_df in df.groupby('residue_name'):
    needed_atoms = expected_atoms[resn]
    needed_atoms.sort()
    # for others, no order needed_atoms.sort()
    for resi,resi_df in resn_df.groupby('residue_number'):
        atoms = resi_df.atom_name.to_list()
        if resi_df.residue_number.iloc[0] == first_res:
            atoms = atoms[:-1] + ["O5'","O6",'OP1','OP2','P']
        if atoms != needed_atoms: # for others shoudl be sorted(atoms)
            print("ff",resi)

# predictions
# checked rna atom roder
for pdb in tqdm(glob(f'{args.input_folder}/rna_R1260TS???.csv')):

    if os.path.isfile(f"{pdb[:-4]}_reordered.csv"):
        print(f'already processed {pdb}')
        continue
    df = pd.read_csv(pdb)
    
    # rename some ambigious atoms
    df = df.replace({'atom_name':{'O1P':'OP1','O2P':'OP2'}})
    
    # confirm number starts correctly
    first_res = df.residue_number.min()
    if first_res == 1:
        df['residue_number'] = df['residue_number']+21
    first_res = df.residue_number.min()
    if first_res != 22:
        print('not 22 start?',pdb)

    # reorder and save
    df = df.sort_values(by=['frame','residue_number','atom_name']) 
    # remove 5' phosphate
    df = df[~((df.residue_number==22)&(df.atom_name.isin(["OP3","O5'",'OP1','OP2','P'])))]
    df.to_csv(f"{pdb[:-4]}_reordered.csv", index=False)
    

    # check content and order
    for (resn,frame),resn_df in df.groupby(['residue_name','frame']):
        needed_atoms = expected_atoms[resn]
        needed_atoms.sort()
        for resi,resi_df in resn_df.groupby('residue_number'):
            atoms = resi_df.atom_name.to_list()
            if resi_df.residue_number.iloc[0] == first_res:
                atoms = atoms[:-1] + ["O5'",'O6','OP1','OP2','P']
            if sorted(atoms) != needed_atoms:
                print("gg",pdb,frame, resi,atoms,needed_atoms)

# no need to check solvent, their identifies are already dealt with in 0
