import numpy as np
import pandas as pd
from glob import glob
from tqdm import tqdm
import os
import argparse
from scipy.spatial.transform import Rotation
from scipy.spatial.distance import cdist

'''
python 2-casp16_water-get_areas_of_interest.py \
    --native_csv 9CBU-rna.csv \
    --input_folder R1260_serial \
    --res 100 \
    --neighborhood_radius 12 \
    --solvent_radius 5 \
    --neighborhood_output neighborhoods \
    --alignment_method all_heavy_atom \
     --num_samples 5
'''

parser = argparse.ArgumentParser(description="Get neighorhoods around each RNA nucleotide.")
parser.add_argument('--native_csv', required=True, help="The native pdb containing rna to align to aka reference.")
parser.add_argument('--input_folder', required=True, help="Folder of the predicted CSV files.")
#parser.add_argument('--group', required=True, type=str, help="The group to process.")
parser.add_argument('--res', required=True, type=int, help="The residue to process.")
parser.add_argument('--neighborhood_radius', required=True, type=float, help="The radius of neighorhood.")
parser.add_argument('--solvent_radius', required=True, type=float, help="The radius of neighorhood.")
parser.add_argument('--neighborhood_output', required=True, help="The fodler to output neighborhoods to.")
parser.add_argument('--alignment_method', required=True, help="The alignment_method, all_heavy_atom, backbone, 3_atom, 5_atom.")
parser.add_argument('--num_samples', required=True, type=int, help="Number of samples, with replacement, to make.")
args = parser.parse_args()

def dist_minRef_to_all(ref_coords, coords):
    #dists = np.linalg.norm(ref_coords[:, np.newaxis, :] - coords[np.newaxis, :, :], axis=-1)
    dists = cdist(ref_coords, coords, metric='euclidean')
    min_dists = dists.min(axis=0)
    return min_dists

def choose_alignment_atoms(rna,alignment_method):
    if alignment_method == 'all_heavy_atom':
        return rna[['x_coord','y_coord','z_coord']].to_numpy().astype(float)
    elif alignment_method == 'backbone':
        backbone_atoms = ['P','OP1','OP2',"O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'"]
        selected_atoms = rna[rna.atom_name.isin(backbone_atoms)]
        return selected_atoms[['x_coord','y_coord','z_coord']].to_numpy().astype(float)
    elif alignment_method == '3_atom':
        atom_3 = [('G','P'),('G',"C4'"),('G','N9'),
                ('A','P'),('A',"C4'"),('A','N9'),
                ('C','P'),('C',"C4'"),('C','N1'),
                ('U','P'),('U',"C4'"),('U','N1'),]
        selected_atoms = rna[rna[['residue_name', 'atom_name']].apply(tuple, axis=1).isin(atom_3)]
        return selected_atoms[['x_coord','y_coord','z_coord']].to_numpy().astype(float)
    elif alignment_method == '5_atom':
        atom_5 = [('G','P'),('G',"C4'"),('G','N9'),('G','C2'),('G','C6'),
                ('A','P'),('A',"C4'"),('A','N9'),('A','C2'),('A','C6'),
                ('C','P'),('C',"C4'"),('C','N1'),('C','C2'),('C','C4'),
                ('U','P'),('U',"C4'"),('U','N1'),('U','C2'),('U','C4')]
        selected_atoms = rna[rna[['residue_name', 'atom_name']].apply(tuple, axis=1).isin(atom_5)]
        return selected_atoms[['x_coord','y_coord','z_coord']].to_numpy().astype(float)
    else:
        print('Alginement method not recongnized.')
        return None


def csv2pdb(df,pdb_f):
    with open(pdb_f,'w') as f:
        for i,atom in df.iterrows():
             atom_id = ' '*(5-len(str((i%100000)))) + str((i%100000))
             if (type(atom.atom_name) is not str) and (type(atom.element_symbol) is not str) and (type(atom.residue_name) is not str):
                 atom_name = "CL  "
                 res_name = " CL"
                 elem = "CL"
             else:
                 if type(atom.atom_name) is not str:
                     atom_name = atom.element_symbol.upper() + ' '*(4-len(atom.element_symbol))
                 else:
                     atom_name = atom.atom_name + ' '*(4-len(atom.atom_name))
                 if type(atom.residue_name) is not str:
                     if type(atom.element_symbol) is not str:
                         res_name = ' '*(3-len(atom.atom_name)) + atom.atom_name.upper()
                     else:
                         res_name = ' '*(3-len(atom.element_symbol.upper())) + atom.element_symbol.upper()
                 else:
                     res_name = ' '*(3-len(atom.residue_name)) + atom.residue_name
                 try:
                     elem = atom.element_symbol + ' '*(2-len(atom.element_symbol))
                 except:
                     elem = f'{atom.atom_name[0]} '
             if atom.residue_number != atom.residue_number: 
                 resnum = '   1'
             else:
                 resnum = ' '*(4-len(str(int(atom.residue_number)))) + str(int(atom.residue_number))
             x = str(round(atom.x_coord,3))
             x += '0'*(3-len(x.split('.')[1]))
             y = str(round(atom.y_coord,3))
             y += '0'*(3-len(y.split('.')[1]))
             z = str(round(atom.z_coord,3))
             z += '0'*(3-len(z.split('.')[1]))
             coords = ' '*(8-len(x)) + x
             coords += ' '*(8-len(y)) + y
             coords += ' '*(8-len(z)) + z
             pdb_line = f'ATOM  {atom_id} {atom_name} {res_name} 0{resnum}    {coords}  1.00  0.00           {elem} \n'
             f.write(pdb_line)

if not os.path.isdir(args.neighborhood_output):
    try: os.makedirs(args.neighborhood_output)
    except FileExistsError: pass
if not os.path.isdir(f'{args.neighborhood_output}/reference'):
    try: os.makedirs(f'{args.neighborhood_output}/reference')
    except FileExistsError: pass
if not os.path.isdir(f'{args.neighborhood_output}/{args.alignment_method}'):
    try: os.makedirs(f'{args.neighborhood_output}/{args.alignment_method}')
    except FileExistsError: pass

# Reference
df = pd.read_csv(args.native_csv)

# get RNA nucleotides in reference neighborhood_radius away
coords = df[['x_coord','y_coord','z_coord']].to_numpy().astype(float)
ref_croods = df[df.residue_number==args.res][['x_coord','y_coord','z_coord']].to_numpy().astype(float)
dists = dist_minRef_to_all(ref_croods,coords)
residues = df[dists<args.neighborhood_radius].residue_number.unique()

# save the neighborhood
if not os.path.isfile(f'{args.neighborhood_output}/reference/{args.res}_9CBU-rna.csv'):
    df[df.residue_number.isin(residues)].to_csv(f'{args.neighborhood_output}/reference/{args.res}_9CBU-rna.csv',index=False)
# get the reference coordinates that will be aligned to
ref_coords = choose_alignment_atoms(df[df.residue_number.isin(residues)],args.alignment_method) #ref[['x_coord','y_coord','z_coord']].to_numpy().astype(float)
ref_center = ref_coords.mean(axis=0)
ref_coords_centered = ref_coords - ref_center
rmsd_constant = np.sqrt(1/ref_coords_centered.shape[0])
csv2pdb(df,f'{args.neighborhood_output}/reference/{args.res}_9CBU-rna.pdb')

results_df = pd.DataFrame(columns=['file', 'frame', 'rmsd','group','residue'])


for f in glob(f'{args.input_folder}/rna_R1260TS*_reordered.csv'):
    # figure out group and initialize its data if needed
    group = f.split('TS')[-1].split('_')[0]

    out_dir = f'{args.neighborhood_output}/{args.alignment_method}/{group}'
    if os.path.isfile(f'{out_dir}/{args.res}_rna.csv'):
        print(f'already processed {group} for res {args.res}')
        continue

    #if group not in sol_dfs:
    #    sol_dfs[group],rna_dfs[group] = [],[]
    sol_dfs = []
    rna_dfs = []
    # get the rna and solvent coords
    rna_df = pd.read_csv(f)
    sol_df = pd.read_csv(f'{"/".join(f.split("/")[:-1])}/sol{f.split("/")[-1][3:-14]}.csv')
    
    # get the solvent that is within solvent_radius of the predicted residue in this frame
    print(group)
    for frame,rna_df_frame in tqdm(rna_df.groupby('frame')):
        sol_df_frame = sol_df[sol_df.frame==frame]
        rna_df_frame_area = rna_df_frame[rna_df_frame.residue_number.isin(residues)].copy()
        sol_ref_coords = rna_df_frame_area[['x_coord','y_coord','z_coord']].to_numpy().astype(float)
        sol_coords = sol_df_frame[['x_coord','y_coord','z_coord']].to_numpy().astype(float)
        dists = dist_minRef_to_all(sol_ref_coords,sol_coords)
        sol_df_frame_area = sol_df_frame[dists<args.solvent_radius].copy()#sol_df.residue_number.isin(sols)]
        

        rna_coords = choose_alignment_atoms(rna_df_frame_area,args.alignment_method)
        all_rna_coords = rna_df_frame_area[['x_coord','y_coord','z_coord']].to_numpy().astype(float)
        sol_coords = sol_df_frame_area[['x_coord','y_coord','z_coord']].to_numpy().astype(float)
        rna_center = rna_coords.mean(axis=0)
        rna_coords_centered = rna_coords - rna_center

        # Check if lengths of vectors are different
        # specifically group 006 missed residue 408 in 10 models for some reason
        if len(ref_coords_centered) != len(rna_coords_centered):
            print(f"Error: Vectors are of different lengths at frame {frame} group {group}")
            continue 
        # align based on rna, only selected atoms
        rot, rssd = Rotation.align_vectors(ref_coords_centered,rna_coords_centered)
        # rssd is the square root of the sum of squares (the weights default to one), rmsd will just differ by sqrt(1/N) where N is number atoms
        rmsd = rssd * rmsd_constant
        
        # align all rna atoms
        rna_aligned = rot.apply(all_rna_coords-rna_center)+ref_center #rna_coords_centered) + ref_center
        # then apply alignment to solvent
        sol_aligned = rot.apply(sol_coords-rna_center) + ref_center
        rna_df_frame_area[['x_coord','y_coord','z_coord']] = rna_aligned
        sol_df_frame_area[['x_coord','y_coord','z_coord']] = sol_aligned


        sol_dfs.append(sol_df_frame_area)
        rna_dfs.append(rna_df_frame_area)
        new_row = pd.DataFrame([{'frame': frame, 'rmsd': rmsd,
            'group': group,'residue': args.res}])
        if results_df.empty:
            results_df = new_row  
        else:
            results_df = pd.concat([results_df, new_row], ignore_index=True)


    #for group in rna_dfs.keys():
    out_dir = f'{args.neighborhood_output}/{args.alignment_method}/{group}'
    if not os.path.isdir(out_dir):
        try: os.makedirs(out_dir)
        except FileExistsError: pass

    all_frames_rna = pd.concat(rna_dfs,ignore_index=True)
    all_frames_sol = pd.concat(sol_dfs,ignore_index=True)
    all_frames_rna.to_csv(f'{out_dir}/{args.res}_rna.csv',index=False)
    all_frames_sol.to_csv(f'{out_dir}/{args.res}_sol.csv',index=False)
    csv2pdb(all_frames_rna,f'{out_dir}/{args.res}_rna.pdb')
    csv2pdb(all_frames_sol,f'{out_dir}/{args.res}_sol.pdb')

    frames = list(all_frames_rna.frame.unique())
    for sample in range(1,args.num_samples+1):
        if not os.path.isdir(f'{out_dir}/sample{sample}'):
            try: os.makedirs(f'{out_dir}/sample{sample}')
            except FileExistsError: pass
        sample_frames = np.random.choice(frames,size=len(frames),replace=True)
        sample_rna_df = pd.concat([all_frames_rna[all_frames_rna.frame==frame] for frame in sample_frames],ignore_index=True)
        sample_sol_df = pd.concat([all_frames_sol[all_frames_sol.frame==frame] for frame in sample_frames],ignore_index=True)
        sample_rna_df.to_csv(f'{out_dir}/sample{sample}/{args.res}_rna.csv',index=False)
        sample_sol_df.to_csv(f'{out_dir}/sample{sample}/{args.res}_sol.csv',index=False)
        csv2pdb(sample_rna_df,f'{out_dir}/sample{sample}/{args.res}_rna.pdb')
        csv2pdb(sample_sol_df,f'{out_dir}/sample{sample}/{args.res}_sol.pdb')

if len(results_df)>0:
    rmsd_file = f'{args.neighborhood_output}/{args.alignment_method}/{args.res}_{args.alignment_method}_rmsd.csv'
    if os.path.isfile(rmsd_file):
        results_df = pd.concat([pd.read_csv(rmsd_file),results_df],ignore_index=True)
    results_df.to_csv(rmsd_file, index=False)
