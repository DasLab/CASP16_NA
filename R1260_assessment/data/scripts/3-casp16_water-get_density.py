import MDAnalysis
from MDAnalysis.analysis import density
from scipy.interpolate import RegularGridInterpolator
from glob import glob
from tqdm import tqdm
import pandas as pd
import numpy as np
import shutil
import mrcfile
import time
import warnings
import argparse
import os
import math
from scipy.spatial.distance import cdist
import gc
import numba

'''
python scripts/3-casp16_water-get_density.py \
    --map scripts/Con2-2.2A_sh.mrc \
    --neighborhoods neighborhoods_10 \
    --group 167 \
    --output_folder avg_maps \
    --alignment_method all_heavy_atom
''' 

parser = argparse.ArgumentParser(description="Align RNA and move solvent along with it.")
parser.add_argument('--map', required=True, type=str, help="The cryoEM map.")
parser.add_argument('--neighborhoods', required=True, help="The neighborhood folder.")
parser.add_argument('--group', required=True, type=str, help="The group to process.")
parser.add_argument('--output_folder', required=True, type=str, help="Box.")
parser.add_argument('--alignment_method', required=True, help="The alignment_method, all_heavy_atom, backbone, 3_atom, 5_atom.")
parser.add_argument('--sample', required=False, type=int, default=None, help="The alignment_method, all_heavy_atom, backbone, 3_atom, 5_atom.")
args = parser.parse_args()

# Ignore specific user warnings
warnings.filterwarnings("ignore", category=UserWarning, message="Box padding .* is not used in user defined grids.")
warnings.filterwarnings("ignore", category=UserWarning, message="Reader has no dt information, set to 1.0 ps")
warnings.filterwarnings("ignore", category=UserWarning, message="Unknown element M found for some atoms.")
warnings.filterwarnings("ignore", category=UserWarning, message="Unknown entry g encountered in formal charge field.")
warnings.filterwarnings("ignore", category=UserWarning, message="Failed to guess the mass for the following atom types:")
warnings.filterwarnings("ignore", category=UserWarning, message="Unknown entry l encountered in formal charge field.")
warnings.filterwarnings("ignore", category=UserWarning, message="No atoms in AtomGroup at input time frame.")
warnings.filterwarnings("ignore", category=UserWarning, message="Unknown entry a encountered in formal charge field.")


@numba.njit(fastmath=True, parallel=True)
def calculate_density_sum_numba(dist2_chunk, pref_chunk, invs2_chunk):
    """
    Calculates the sum of gaussians using Numba for maximum speed.
    
    Args:
        dist2_chunk (array): Shape (N_atoms, N_voxels)
        pref_chunk (array): Shape (N_atoms, 5)
        invs2_chunk (array): Shape (N_atoms, 5)
        
    Returns:
        array: Shape (N_voxels,)
    """
    n_atoms, n_voxels = dist2_chunk.shape
    n_gaussians = pref_chunk.shape[1]
    
    # Initialize the output array for the voxels in this chunk
    map_sum_chunk = np.zeros(n_voxels, dtype=np.float32)

    # Numba will parallelize this outer loop automatically
    for j in numba.prange(n_voxels):
        voxel_density = 0.0
        for i in range(n_atoms):
            dist2 = dist2_chunk[i, j]
            for k in range(n_gaussians):
                exponent = -0.5 * invs2_chunk[i, k] * dist2
                voxel_density += pref_chunk[i, k] * np.exp(exponent)
        map_sum_chunk[j] = voxel_density
        
    return map_sum_chunk


def create_min_distance_grid(xs, ys, zs, res_coords, dist_cutoff):

    x_min, x_max = res_coords[:, 0].min() - dist_cutoff, res_coords[:, 0].max() + dist_cutoff
    y_min, y_max = res_coords[:, 1].min() - dist_cutoff, res_coords[:, 1].max() + dist_cutoff
    z_min, z_max = res_coords[:, 2].min() - dist_cutoff, res_coords[:, 2].max() + dist_cutoff

    xs_indices = (xs >= x_min) & (xs <= x_max)
    ys_indices = (ys >= y_min) & (ys <= y_max)
    zs_indices = (zs >= z_min) & (zs <= z_max)

    valid_xs, valid_ys, valid_zs = xs[xs_indices], ys[ys_indices], zs[zs_indices]

    X, Y, Z = np.meshgrid(valid_xs, valid_ys, valid_zs, indexing='ij')
    grid_points = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=-1)

    dists = cdist(grid_points, res_coords, metric='euclidean')

    min_dists = np.min(dists, axis=1).reshape(X.shape)

    min_distance_grid = np.full((len(xs), len(ys), len(zs)), np.nan)
    
    mask = min_dists < dist_cutoff

    valid_region = np.ix_(xs_indices, ys_indices, zs_indices)
    min_distance_grid[valid_region] = min_dists  # Assign all distances first
    min_distance_grid[valid_region] = np.where(mask, min_dists, np.nan)

    return min_distance_grid


def create_density_grid(xs, ys, zs, get_density, min_distance_grid):
    density_grid = np.full(min_distance_grid.shape, np.nan)
    valid_indices = ~np.isnan(min_distance_grid)

    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
    grid_points = np.stack([X[valid_indices], Y[valid_indices], Z[valid_indices]], axis=-1)

    densities = get_density(grid_points)
    density_grid[valid_indices] = densities

    return density_grid


def calculate_weighted_average_density(density_grids, min_distance_grids):
    sum_weighted_densities = np.zeros(density_grids[0].shape)
    sum_weights = np.zeros(density_grids[0].shape)
    
    for density_grid, min_distance_grid in zip(density_grids, min_distance_grids):
        valid_indices = ~np.isnan(min_distance_grid)
        weights = np.zeros(min_distance_grid.shape)
        weights[valid_indices] = 1 / min_distance_grid[valid_indices]

        sum_weighted_densities += np.nan_to_num(density_grid) * weights
        sum_weights += weights
    
    weighted_avg_density = np.divide(sum_weighted_densities, sum_weights, where=sum_weights != 0)
    return weighted_avg_density

def get_density_interpolator_dx(grid):
    xs,ys,zs = grid.midpoints
    mapvalues = grid.grid

    get_mapvalue = RegularGridInterpolator((xs,ys,zs),mapvalues, method='linear')
    return get_mapvalue

def get_density_for_pdb(rna_out_name,grid_center,xdim,ydim,zdim,
                        selection='all',
                       step_size=1,pixel_size=0.82,):
    rna = MDAnalysis.Universe(rna_out_name).select_atoms(selection)
    dens = density.DensityAnalysis(rna, delta=pixel_size, 
        gridcenter=grid_center, xdim=xdim, ydim=ydim, zdim=zdim) 
    dens.run(step=step_size, verbose=False)
    dens.results.density.convert_density('water')
    get_density = get_density_interpolator_dx(dens.results.density)
    return get_density, dens.results.density, rna

### SIGMA
B_={}
B_["C"]=np.array([0.2465, 1.7100, 6.4094, 18.6113, 50.2523])
B_["O"]=np.array([0.2067, 1.3815, 4.6943, 12.7105, 32.4726])
B_["N"]=np.array([0.2451, 1.7481, 6.1925, 17.3894, 48.1431])
B_["S"]=np.array([0.2681, 1.6711, 7.0267, 19.5377, 50.3888])
B_["P"]=np.array([0.2908, 1.8740, 8.5176, 24.3434, 63.2996])
B_["F"]=np.array([0.2057, 1.3439, 4.2788, 11.3932, 28.7881])
B_["MG"]=np.array([0.3278,2.2720, 10.9241, 39.2898, 101.9748])
B_["K"]=np.array([0.3703, 3.3874, 13.1029, 68.9592, 194.4329])
B_["NA"]=np.array([0.3334, 2.3446, 10.0830, 48.3037, 138.2700])
B_["CL"]=np.array([0.2468, 1.5242, 6.1537, 16.6687, 42.3086])
### WEIGHT
A_={}
A_["C"]=np.array([0.0893, 0.2563, 0.7570, 1.0487, 0.3575])
A_["O"]=np.array([0.0974, 0.2921, 0.6910, 0.6990, 0.2039])
A_["N"]=np.array([0.1022, 0.3219, 0.7982, 0.8197, 0.1715])
A_["S"]=np.array([0.2497, 0.5628, 1.3899, 2.1865, 0.7715])
A_["P"]=np.array([0.2548, 0.6106, 1.4541, 2.3204, 0.8477])
A_["F"]=np.array([0.1083, 0.3175, 0.6487, 0.5846, 0.1421])
A_["MG"]=np.array([0.2314, 0.6866, 0.9677, 2.1882,1.1339])
A_["K"]=np.array([0.4115, 1.4031, 2.2784, 2.6742, 2.2162])
A_["NA"]=np.array([0.2142, 0.6853, 0.7692, 1.6589, 1.4482])
A_["CL"]=np.array([0.2443, 0.5397, 1.3919, 2.0197, 0.6621])

def get_fmod_param(atoms):
    pref=[]; invs2=[]
    # cycle on atoms
    for at in atoms:
        # get atom type
        atype = at.type

        # RCK hacky way to deal with problems with PDB
        if atype == "M":
            atype = "MG"
        elif atype == "N" and at.resname == 'N':
            atype = "NA"

        # add to pref list
        pref.append(A_[atype] * pow(math.pi / ( B_[atype] + at.tempfactor/4.0 ), 1.5))
        # add to invs2 list
        invs2.append(2.0 * math.pi * math.pi / ( B_[atype] + at.tempfactor/4.0 ))
    # convert to tensor and copy to device [natoms, 5]
    pref  = np.array(pref)
    invs2 = np.array(invs2)
    return pref, invs2

def get_density_grid_scatter(atoms, density_grid_other, xs, ys, zs, cutoff=4):
    pref, invs2 = get_fmod_param(atoms)  # Natoms, 5

    density_grid = np.full(density_grid_other.shape, np.nan)
    valid_indices = ~np.isnan(density_grid_other)

    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
    grid_points = np.stack([X[valid_indices], Y[valid_indices], Z[valid_indices]], axis=-1)
    
    pos_g = atoms.positions
    num_atoms = pos_g.shape[0]
    num_voxels = grid_points.shape[0]
    
    map_sum = np.zeros(num_voxels)
    

    # Spatial sorting using linear
    round_dist = 2
    rounded_pos_g = (np.round(pos_g/round_dist)*round_dist).astype(int)
    sorted_indices = np.lexsort((rounded_pos_g[:, 2], rounded_pos_g[:, 1], rounded_pos_g[:, 0]))
    sorted_positions = pos_g[sorted_indices]
    pref = pref[sorted_indices]
    invs2 = invs2[sorted_indices]

    atom_chunk_size = max(1,len(pos_g)//1400) # TODO option means nothing now
    for atom_start in range(0, num_atoms, atom_chunk_size):
        atom_end = min(atom_start + atom_chunk_size, num_atoms)
        atom_chunk = sorted_positions[atom_start:atom_end]

        x_min = atom_chunk[:, 0].min() - cutoff
        x_max = atom_chunk[:, 0].max() + cutoff
        y_min = atom_chunk[:, 1].min() - cutoff
        y_max = atom_chunk[:, 1].max() + cutoff
        z_min = atom_chunk[:, 2].min() - cutoff
        z_max = atom_chunk[:, 2].max() + cutoff

        within_cutoff_mask = ((grid_points[:, 0] >= x_min) & (grid_points[:, 0] <= x_max) &
            (grid_points[:, 1] >= y_min) & (grid_points[:, 1] <= y_max) &
            (grid_points[:, 2] >= z_min) & (grid_points[:, 2] <= z_max))
        chunk_indices = np.where(within_cutoff_mask)[0]
        valid_vox_chunk = grid_points[chunk_indices]

        dist2_chunk = cdist(atom_chunk, valid_vox_chunk, metric='sqeuclidean')
        # cutoff things more
        within_cutoff_dist_mask = (dist2_chunk <= cutoff ** 2).any(axis=0)
        chunk_indices = chunk_indices[within_cutoff_dist_mask]
        valid_vox_chunk = valid_vox_chunk[within_cutoff_dist_mask]
        dist2_chunk = dist2_chunk[:, within_cutoff_dist_mask]

        map_sum_chunk = calculate_density_sum_numba(dist2_chunk, pref[atom_start:atom_end], invs2[atom_start:atom_end])

        map_sum[chunk_indices] += map_sum_chunk

    density_grid[valid_indices] = map_sum

    return density_grid




def submap_average(pdb_files,selection,mrc_to_map_to,output,dist_cutoff,box_radius,ref_location):
    # for each submap
    submaps = []
    coords_list = []
    sum_weighted_densities = None
    pdbs = glob(pdb_files)

    if len(pdbs) == 0:
        print(f'No pdbs {pdb_files}')
        return None

    with mrcfile.open(mrc_to_map_to, mode='r') as mrc:
        xs = (np.arange(mrc.header.nxstart,mrc.header.mx)*mrc.voxel_size.x)+mrc.header.origin.x
        ys = (np.arange(mrc.header.nystart,mrc.header.my)*mrc.voxel_size.y)+mrc.header.origin.y
        zs = (np.arange(mrc.header.nzstart,mrc.header.mz)*mrc.voxel_size.z)+mrc.header.origin.z
        zero_shape = mrc.data.shape
    for pdb_file in tqdm(pdbs,desc='loading in densities'):
        # get model
        res_num = pdb_file.split('/')[-1].split('_')[0]
        ref_csv_file = f'{ref_location}/{res_num}_9CBU-rna.csv'
        ref_pdb_file = f'{ref_location}/{res_num}_9CBU-rna.pdb'
        # extract nucleotide coordintes
        # this is only the res_num of interest
        pdb = pd.read_csv(ref_csv_file)
        pdb = pdb[pdb.residue_number == int(res_num)]
        res_coords = pdb[['x_coord','y_coord','z_coord']].to_numpy().astype(float)

        ref = MDAnalysis.Universe(ref_pdb_file)
        x,y,z = ref.select_atoms('all').positions.T
        xdim = x.max()-x.min()+(2*box_radius)
        ydim = y.max()-y.min()+(2*box_radius)
        zdim = z.max()-z.min()+(2*box_radius)
        grid_center = (x.min()+((x.max()-x.min())/2), y.min()+((y.max()-y.min())/2), z.min()+((z.max()-z.min())/2))
        # get interpolation object
        if os.path.getsize(pdb_file) == 0:
            continue # no atoms, empty pdb
        get_density, grid, atoms = get_density_for_pdb(pdb_file,grid_center,xdim,ydim,zdim,selection)

        min_distance_grid = create_min_distance_grid(xs, ys, zs, res_coords, dist_cutoff)
        density_grid = create_density_grid(xs, ys, zs, get_density, min_distance_grid)
        
        
        # get_scatter, scatter_dens = get_density_for_pdb_scatter(atoms,density_grid)
        scatter_grid = get_density_grid_scatter(atoms,density_grid,xs, ys, zs)
        
        if sum_weighted_densities is None:
            sum_weighted_densities = np.zeros(density_grid.shape)
            sum_weights = np.zeros(density_grid.shape)
            sum_weighted_scatter = np.zeros(density_grid.shape)
        valid_indices = ~np.isnan(min_distance_grid)
        weights = np.zeros(min_distance_grid.shape)
        weights[valid_indices] = 1 / min_distance_grid[valid_indices]

        sum_weighted_densities += np.nan_to_num(density_grid) * weights
        sum_weighted_scatter += np.nan_to_num(scatter_grid) * weights
        sum_weights += weights

    if sum_weighted_densities is None:
        sum_weighted_densities = np.zeros(zero_shape)
        sum_weights = np.zeros(zero_shape)
        sum_weighted_scatter = np.zeros(zero_shape)
    weighted_avg_density = np.divide(sum_weighted_densities, sum_weights, where=sum_weights != 0)
    weighted_avg_scatter = np.divide(sum_weighted_scatter, sum_weights, where=sum_weights != 0)

    shutil.copy(mrc_to_map_to,f'{output}_density.mrc')
    shutil.copy(mrc_to_map_to,f'{output}_scatter.mrc')

    # get 2.2A coordinate frame
    with mrcfile.open(f'{output}_density.mrc', mode='r+') as mrc_dens, mrcfile.open(f'{output}_scatter.mrc', mode='r+') as mrc_scat :
        # Ensure appropriate axes are in 1-to-1 correspondance with MRC format
        if mrc_dens.header.mapc == 1:
            weighted_avg_density = np.swapaxes(weighted_avg_density, 0, 2)
            mrc_dens.header.mapc = 0
            weighted_avg_scatter = np.swapaxes(weighted_avg_scatter, 0, 2)
            mrc_scat.header.mapc = 0
        
        mrc_dens.data[:] = weighted_avg_density
        mrc_scat.data[:] = weighted_avg_scatter


if args.sample: 
    output = f'{args.output_folder}/{args.neighborhoods}_{args.alignment_method}/sample{args.sample}'
    folder = f'{args.neighborhoods}/{args.alignment_method}/{args.group}/sample{args.sample}'
else: 
    folder = f'{args.neighborhoods}/{args.alignment_method}/{args.group}'
    output = f'{args.output_folder}/{args.neighborhoods}_{args.alignment_method}'
if not os.path.isdir(output):
    os.makedirs(output)

print('processing RNA')
if not os.path.isfile(f'{output}/{args.group}_rna_scatter.mrc'):
    submap_average(pdb_files=f'{folder}/*_rna.pdb',
        selection='all',
        mrc_to_map_to=args.map,
        output=f'{output}/{args.group}_rna',
        dist_cutoff=20,
        box_radius=30,
        ref_location=f'{args.neighborhoods}/reference/')

sol_selections = {'wat':'resname HOH or resname WAT or resname SOL',
                  'mg':'resname MG or resname nMg or resname mMg',
                  'na':'resname NA or resname Na+',
                  'cl':'resname CL or resname Cl- or resname UNL',
                  'k':'resname K'}

for mol,selection in sol_selections.items():
    print('processing',mol)
    if not os.path.isfile(f'{output}/{args.group}_{mol}_scatter.mrc'):
        submap_average(pdb_files=f'{folder}/*_sol.pdb',
            selection=selection,
            mrc_to_map_to=args.map,
            output=f'{output}/{args.group}_{mol}',
            dist_cutoff=20,
            box_radius=30,
            ref_location=f'{args.neighborhoods}/reference/')
