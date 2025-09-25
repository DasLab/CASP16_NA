# imports
from os import system
from glob import glob
from io import StringIO
import numpy as np
from tqdm import tqdm
import pandas as pd


###############################################################################
# Plotting 
###############################################################################
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from matplotlib.gridspec import GridSpec
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker 
font_path = font_manager.findfont(font_manager.FontProperties(family='Arial'))
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = [font_path]
plt.rcParams['svg.fonttype'] = 'none'

###############################################################################
# Shared varibales
###############################################################################

# Target lists
###############################################################################
ALL_MONOMER = ['R0250', 'R0251', 'R0252', 'R0253', 'R0253v1', 'R0253v2', 'R0254', 'R0281','R0283', 'R0285', 'R0290', 
              'R1203', 'R1205', 'R1209', 'R1211', 'R1212', 'R1221s2',
              'R1221s3', 'R1224s2', 'R1224s3', 'R1241', 'R1242', 'R1248', 'R1250', 'R1251',
              'R1252', 'R1253v1', 'R1253v2', 'R1254', 'R1255', 'R1256', 'R1261', 'R1262',
              'R1263', 'R1264', 'R1271', 'R1281', 'R1283v1', 'R1283v2','R1283v3',
               'R1285', 'R1286', 'R1288', 'R1289',
              'R1290', 'R1291', 'R1293', 'R1296', 'D1273']
MONOMER_TARGETS_TO_SCORE = ALL_MONOMER.copy()

# for monomer assessement we are not scoring round 0
DO_NOT_SCORE_MONOMER = ['R0250', 'R0251', 'R0252', 'R0253', 'R0253v1', 'R0253v2', 'R0254', 'R0281', 'R0283', 'R0285', 'R0290'] #, 'R1262']
[MONOMER_TARGETS_TO_SCORE.remove(t) for t in DO_NOT_SCORE_MONOMER]
MONOMER_TARGETS_TO_SCORETEMP = MONOMER_TARGETS_TO_SCORE.copy()
[MONOMER_TARGETS_TO_SCORETEMP.remove(t) for t in ['R1283v2','R1283v3']]
# these are all there similiar (same seq) to only count one score to final score
TARGETS_CHOOSE_BEST = [['R1253v1', 'R1253v2'],['R1263', 'R1264'],['R1261', 'R1262'],
                                ['R1283v1', 'R1283v2','R1283v3'],['R1241', 'R1291'],['R1221s3', 'R1224s3'],
                                ['M1228v1', 'M1228v2'],
                 ['M1239v1', 'M1239v2']]


MULTIMER_TARGETS = ['R1250', 'R1251', 'R1252', 'R1253v1', 'R1253v2',
       'R1254', 'R1281', 'R1283v2', 'R1283v3', 'R1285', 'R1290']
#MULTIMER_TARGETS_CHOOSE_BEST = [['R1253v1', 'R1253v2'],['R1283v2','R1283v3']]


HYBRID_TARGETS = ['M1209', 'M1211', 'M1212', 'M1216', 
                 'M1221', 'M1224', 'M1228v1', 'M1228v2',
                 'M1239v1', 'M1239v2', 'M1268', 'M1271', 
                  'M1276', 'M1282', 'M1287', 'M1293', 
                  'M1296', 'M1297'] 
#HYBRID_TARGETS_CHOOSE_BEST = [['M1228v1', 'M1228v2'],
#                 ['M1239v1', 'M1239v2']]

LIGAND_TARGETS = ['D1273','R1261','R1263','R1264',  'R1288']
#LIGAND_TARGETS_CHOOSE_BEST = [['R1263', 'R1264']]
# R1262 excluded because incorrect ligand

# Metrics data
###############################################################################

# if applicable
METRIC_PASS_SCORE = {'tm_score':0.45,
                    'gdt_ts':0.45,
                    'lddt':0.75,
                    'lddt_no_checks':0.75,
                    'lddt_pli':0.5,
                    'lddt_lp':0.5,
                    'ICS(F1)':0.5,'IPS':0.5,'ilDDT':0.5,
                    'NA-NA_ICS(F1)':0.5,'NA-NA_IPS':0.5,'NA-NA_ilDDT':0.5,
                    'NA-protein_ICS(F1)':0.5,'NA-protein_IPS':0.5,'NA-protein_ilDDT':0.5,
                    'f1_bp_list':0.5, 'recall_bp_list':0.5, 
               'precision_bp_list':0.5,'recall_singlets_as_bp':0.5,
                  'f1_singlet':0.5, 'f1_crossed':0.5,'F1':0.5}
METRIC_RANGES = {'global_rmsd':(0,200),
                'clashscore':(0,1000),
                'clashes_number':(0,1000),
                    'rmsd':(0,100),
                'bb_rmsd':(0,20)}

ALL_SCORES = {'Z_CASP15':{"tm_score": 1 / 3, "gdt_ts": 1 / 3, "inf_all":1/8, "lddt_no_checks": 1/8, "clashscore": 1/12},
              'Z_CASP16_t33_g33_l33':{"tm_score": 1 / 3, "gdt_ts": 1 / 3, "lddt": 1/3},
              'Z_CASP16_t30_g30_l40':{"tm_score": 0.3, "gdt_ts": 0.3, "lddt": 0.4},
              'Z_CASP16_t25_g25_l50':{"tm_score": 1 / 4, "gdt_ts": 1 / 4, "lddt": 1/2},
              'Z_CASP16_t20_g20_l60':{"tm_score": 0.2, "gdt_ts": 0.2, "lddt": 0.6},
              'Z_topo':{"tm_score": 1 / 2, "gdt_ts": 1 / 2},
             'Z_rmsd_inf':{"inf_all": 1 / 2, "global_rmsd": 1 / 2}}

# put in the scatters of raw scores --> size issues highlight
# explain lddt with stereo
# plot mean lddt gdt etc (without z) --> ranking based on this
# summary how many past threshold compared to last year
# TM v homolog particulary 

METRICS = {'lddt':'max','tm_score':'max','gdt_ts':'max','inf_all':'max',
           'TM_align':'max',
           'inf_stack':'max','inf_wc':'max','inf_nwc':'max','lddt_no_checks':'max',
          'global_rmsd':'min','clashscore':'min','clashes_number':'min'}
METRICS_LABELS = {'lddt':'lDDT\nsterics checked','tm_score':'TM-score','tm_align':'TM-align','gdt_ts':'GDT-TS','inf_all':'INF',
              'inf_stack':'inf_stack','inf_wc':'inf_wc','inf_nwc':'inf_nwc',
               'lddt_no_checks':'lDDT\nsterics not checked',
          'global_rmsd':'RMSD','clashscore':'clashscore','clashes_number':'clashes'
}
MULTIMER_METRICS = {'lddt':'max','tm_score':'max','gdt_ts':'max',
           'TM_align':'max',
          'global_rmsd':'min','clashes_number':'min',
                   'ICS(F1)':'max', 
                    'Prec.Iface':'max', 
                    'Recal.Iface':'max', 
                    'IPS':'max', 
                    'QSglob':'max', 
                    'QSbest':'max', 
                    'ilDDT':'max', 
                    'GlobDockQ':'max', 
                    'BestDockQ':'max'}
HYBRID_METRICS = MULTIMER_METRICS.copy()
for metric in ['IPS','ICS(F1)','ilDDT']:
    HYBRID_METRICS[metric] = MULTIMER_METRICS[metric]
    for int_type in ['NA-protein','NA-NA','protein-protein','all',"NA-containing"]:
        for method in ['naive_','','log10_']:
            HYBRID_METRICS[f'{int_type}_{method}{metric}'] = MULTIMER_METRICS[metric]

LIGAND_METRICS = {'lddt_pli':'max', 'lddt_lp':'max','rmsd':'min', 'bb_rmsd':'min'}

def _max_with_nan_per_column(x):
    return x.max() if not x.isna().any() else np.nan

def _min_with_nan_per_column(x):
    return x.min() if not x.isna().any() else np.nan
SS_METRICS = {'f1_bp_list':_max_with_nan_per_column, 'recall_bp_list':_max_with_nan_per_column, 
               'precision_bp_list':_max_with_nan_per_column,'recall_singlets_as_bp':_max_with_nan_per_column,
                  'f1_singlet':_max_with_nan_per_column, 'f1_crossed':_max_with_nan_per_column, }
SS_METRICS2 = {'f1_bp_list':_max_with_nan_per_column, 'recall_bp_list':_max_with_nan_per_column, 
               'precision_bp_list':_max_with_nan_per_column }
MOTIF_METRICS = {'F1':_max_with_nan_per_column}

MULTIMER_SCORES = {'Z_interface':{"ICS(F1)": 1 / 3, "IPS": 1 / 3, "ilDDT":1/3}, 
              'Z_topo':{"tm_score": 1 / 2, "gdt_ts": 1 / 2},
                   'Z_monomer':{"tm_score": 0.3, "gdt_ts": 0.3, "lddt": 0.4},
              'Z_0.6_0.4':{"tm_score": 0.6*0.3, "gdt_ts": 0.6*0.3, "lddt": 0.6*0.4,
                           "ICS(F1)": 0.4*1 / 3, "IPS": 0.4*1 / 3, "ilDDT":0.4*1/3},
                  'Z_0.5_0.5':{"tm_score": 0.5*0.3, "gdt_ts": 0.5*0.3, "lddt": 0.5*0.4,
                           "ICS(F1)": 0.5*1 / 3, "IPS": 0.5*1 / 3, "ilDDT":0.5*1/3},
                  'Z_0.4_0.6':{"tm_score": 0.4*0.3, "gdt_ts": 0.4*0.3, "lddt": 0.4*0.4,
                           "ICS(F1)": 0.6*1 / 3, "IPS": 0.6*1 / 3, "ilDDT":0.6*1/3},
                  'Z_0.3_0.7':{"tm_score": 0.3*0.3, "gdt_ts": 0.3*0.3, "lddt": 0.3*0.4,
                           "ICS(F1)": 0.7*1 / 3, "IPS": 0.7*1 / 3, "ilDDT":0.7*1/3},}
HYBRID_SCORES = MULTIMER_SCORES.copy()
for x,y in [(0.6,0.4),(0.5,0.5),(0.4,0.6),(0.3,0.7)]:
    for method in ['naive_','','log10_']:
        HYBRID_SCORES[f'Z_{x}_{y}_NAeven{method}'] = {"tm_score": x*0.3, "gdt_ts": x*0.3, "lddt": x*0.4,
                           f"NA-NA_{method}ICS(F1)": 0.5*y*1 / 3, f"NA-NA_{method}IPS": 0.5*y*1 / 3, f"NA-NA_{method}ilDDT":0.5*y*1/3,
                           f"NA-protein_{method}ICS(F1)": 0.5*y*1 / 3, f"NA-protein_{method}IPS": 0.5*y*1 / 3, f"NA-protein_{method}ilDDT":0.5*y*1/3,}
        HYBRID_SCORES[f'Z_{x}_{y}_NA{method}'] = {"tm_score": x*0.3, "gdt_ts": x*0.3, "lddt": x*0.4,
                           f"NA-containing_{method}ICS(F1)": y*1 / 3, f"NA-containing_{method}IPS": y*1 / 3, f"NA-containing_{method}ilDDT": y*1/3,}



LIGAND_SCORES = {'Z_0.3_0.7':{'lddt_pli':0.3, 'lddt_lp':0.7},
                'Z_0.4_0.6':{'lddt_pli':0.4, 'lddt_lp':0.6},
                'Z_0.5_0.5':{'lddt_pli':0.5, 'lddt_lp':0.5},
                'Z_0.6_0.4':{'lddt_pli':0.6, 'lddt_lp':0.4},
                'Z_0.7_0.3':{'lddt_pli':0.7, 'lddt_lp':0.3}}

MOTIF_TYPES = {'LOOP_E_SUBMOTIF':'Loop-E Submotif', 'INTERCALATED_T_LOOP':'Intercalated T-Loop', 'GNRA_TETRALOOP':'GNRA Tetraloop', 
               'T_LOOP':'T-loop', 'TL_RECEPTOR':'Tetraloop-receptor', 'U_TURN':'U-turn', 'GA_MINOR':'GA-minor', 'UA_HANDLE':'UA-handle', 
               'Z_TURN':'Z-turn', 'PLATFORM':'Platform', 'BULGED_G':'Bulged-G', 'A_MINOR':'A-minor', 'TANDEM_GA_SHEARED':'Tandem-GA-Sheared',
               'INTERCALATED_NT':'Intercalated nucleotide','DOCKED_A':'Docked A'}

# Target and group information
###############################################################################

PARTICIPATION_RATE = 0.6


unresolved_df = pd.read_csv('../raw_scores/unresolved_residues.csv')
UNRESOLVED_RESIDUES = {}
for i,row in unresolved_df.iterrows():
    if pd.isna(row.Unresolved):
        unres = []
    else:
        unres = row.Unresolved.split('-')
    UNRESOLVED_RESIDUES[row.Target] = {'start':row.Start,
                                      'end':row.End,
                                      'unresolved':unres}
UNRESOLVED_RESIDUES['R1203'] = UNRESOLVED_RESIDUES['R1203v1']
UNRESOLVED_RESIDUES['R1255'] = UNRESOLVED_RESIDUES['R1255_1']
UNRESOLVED_RESIDUES['R1256'] = UNRESOLVED_RESIDUES['R1256_1']
UNRESOLVED_RESIDUES['R1281'] = UNRESOLVED_RESIDUES['R1281_fixed0']

#unresolved in bp_list [[85, 100]]



GROUP_INFO = pd.read_csv('../raw_scores/CASP16_groups.csv')
GROUP_INFO['gr_code'] = GROUP_INFO['Group number'].apply(lambda x: '0'*(3-len(str(x)))+str(x))
GROUP_INFO['final_name'] = GROUP_INFO.apply(lambda x: x['Group Name'] + ' - ' + x['gr_code'], axis=1)
GR_CODE_TO_NAME = GROUP_INFO.set_index('gr_code').final_name.to_dict()
GROUP_DICT = GROUP_INFO.set_index('gr_code')['Group Name'].to_dict()
GROUP_SERVER = GROUP_INFO.set_index('final_name')[['Type']].rename(columns={'Type':'Server'})
GR_CODE_TO_NAME_SERVER = GR_CODE_TO_NAME.copy()
for i,name in GROUP_SERVER[GROUP_SERVER.Server=='Server'].iterrows():
    GR_CODE_TO_NAME_SERVER[i[-3:]] = i.rsplit('-',1)[0] + '- S' + i.rsplit('-',1)[1]

TEMPLATE_INFO = pd.read_csv('../raw_scores/historical_overview.csv')
CASP16_TEMPLATE_INFO = TEMPLATE_INFO[TEMPLATE_INFO.Competition == 'CASP16']
TARGET_INFO = pd.read_csv('../raw_scores/casp16_targets.csv')
CASP16_LENGTH = TARGET_INFO.set_index('Target')['Length RNA (nt)'].to_dict()
CASP16_TEMPLATE_DICT = CASP16_TEMPLATE_INFO.set_index('Target')['Best Template TM align'].to_dict()
CASP16_TEMPLATE_DICT['R1283v2'] = CASP16_TEMPLATE_DICT['R1283v1']
CASP16_TEMPLATE_DICT['R1283v3'] = CASP16_TEMPLATE_DICT['R1283v1']
CASP16_TEMPLATE_DICT['R1253v1'] = CASP16_TEMPLATE_DICT['R1253v2']

TEMPLATE_BOUNDARY = 0.45
TEMPLATE_TARGETS = CASP16_TEMPLATE_INFO[CASP16_TEMPLATE_INFO['Best Template TM align']>TEMPLATE_BOUNDARY].Target.to_list()
NO_TEMPLATE_TARGETS = CASP16_TEMPLATE_INFO[CASP16_TEMPLATE_INFO['Best Template TM align']<=TEMPLATE_BOUNDARY].Target.to_list()
NO_TEMPLATE_TARGETS.extend(['R1253v1', 'R1283v2', 'R1283v3'])

NEFF_BOUNDARY = 130
MSA_INFO = pd.read_csv('../raw_scores/msa/all_neff.csv')
CASP16_NEFF = MSA_INFO.set_index('Target')['0.8'].to_dict()
MSA_TARGETS = MSA_INFO[MSA_INFO['0.8']>NEFF_BOUNDARY].Target.to_list()
NO_MSA_TARGETS = MSA_INFO[MSA_INFO['0.8']<=NEFF_BOUNDARY].Target.to_list()
MSA_TARGETS.extend(['R1253v1','R1253v2','R1283v1','R1283v2', 'R1283v3'])
CASP16_NEFF['R1253v1'] = CASP16_NEFF['R1253']
CASP16_NEFF['R1253v2'] = CASP16_NEFF['R1253']
CASP16_NEFF['R1283v1'] = CASP16_NEFF['R1283']


GROUP_HIGHLIGHTS = {GR_CODE_TO_NAME_SERVER['481']:'magenta',
                   GR_CODE_TO_NAME_SERVER['052']:'red',
                   GR_CODE_TO_NAME_SERVER['183']:'orange',
                   GR_CODE_TO_NAME_SERVER['338']:'green',
                   GR_CODE_TO_NAME_SERVER['294']:'cyan',
                   GR_CODE_TO_NAME_SERVER['286']:'blue',
                   GR_CODE_TO_NAME_SERVER['063']:'purple'}

###############################################################################
# code
###############################################################################


# For motif anlalysis
###############################################################################
def is_resolved(residue,target):
    criteria1 = residue >= UNRESOLVED_RESIDUES[target]['start']
    criteria2 = residue <= UNRESOLVED_RESIDUES[target]['end']
    criteria3 = not str(residue) in UNRESOLVED_RESIDUES[target]['unresolved']
    return criteria1 and criteria2 and criteria3





###############################################################################
# secondary structure code
# new parsing from the raw dssr output, of course the ss will be from above.
def parse_base_pair_dssr(row):
    nt1 = int(row.nt1.split('.')[1][1:])
    nt2 = int(row.nt2.split('.')[1][1:])
    chain1 = row.nt1.split('.')[0]
    chain2 = row.nt2.split('.')[0]
    n1 = row.nt1.split('.')[1][0].lower()
    n2 = row.nt2.split('.')[1][0].lower()
    return [f'{chain1}:{n1}{nt1}',f'{chain2}:{n2}{nt2}']

def get_residue_format_from_index(ind,seq,chain=0):
    if isinstance(ind, list):
        return [get_residue_format_from_index(item, seq, chain) for item in ind]
    else:
        return f'{chain}:{seq[ind].lower()}{ind+1}'

    
def convert_dotbracket_to_bp_list(s, allow_pseudoknots=False,silent_values=[]):
    bp_list = []
    lower_alphabet = [chr(lower) for lower in range(97, 123)]
    upper_alphabet = [chr(upper) for upper in range(65, 91)]

    if allow_pseudoknots:
        openDelimiterMap = {
            '(': [],
            '[': [],
            '{': [],
            '<': [],
        }

        # Both the left and right delimiter of each pair of characters point to the same
        # array, so that when we encounter a left character we push an opening-half pair index,
        # and when we encounter a right character we can pop off the most recent one.
        closeDelimiterMap = {
            ')': openDelimiterMap['('],
            ']': openDelimiterMap['['],
            '}': openDelimiterMap['{'],
            '>': openDelimiterMap['<'],
        }

        flipped = None
        for (i, char) in enumerate(s):
            if char == '.':
                # Unpaired base
                continue
            elif (char in openDelimiterMap.keys()):
                # This is an opening delimiter which we've already registered
                openDelimiterMap[char].append(i)
            elif (char in closeDelimiterMap.keys()):
                # This is a closing delimiter which we've already registered               
                partner = closeDelimiterMap[char].pop() if closeDelimiterMap[char] else None
                if (partner == None):
                    raise Exception(f"Unbalanced parenthesis notation: found closing character '{char}'");
                bp_list.append([partner, i])
            elif (
                ((flipped == None or flipped == False) and char.islower())
                or ((flipped == None or flipped == True) and char.isupper())
            ):
                # For performance, we don't initialize our pair stacks with alpha characters.
                # Also, we don't know whether lower/upper case characters are considered left or right
                # delimiters until we encounter one for the first time. Whichever we see first
                # we treat as a left delimiter, then whenever we encounter a new letter, it has
                # to follow the same standard.
                if (flipped == None):
                    flipped = char.isupper()

                pairStack = []
                openDelimiterMap[char] = pairStack 
                closeDelimiterMap[char.lower() if flipped else char.upper()] = pairStack
                pairStack.append(i)
            elif (char.isalpha()):
                # We haven't encountered this character yet, but the case that showed up
                # we have designated as a right delimiter
                raise Exception(f"Unbalanced parenthesis notation: found closing character '{char}'")
            elif char not in silent_values:
                print(f"WARNING: characters in structure, '{char}' ignored!")
                continue
        for (char, pairStack) in openDelimiterMap.items():
            if (len(pairStack)): raise Exception(f"Unbalanced parenthesis notation: found unclosed pair for character '{char}'");

    bp_list = sorted(bp_list, key=lambda x: x[0])
    return bp_list


def convert_bp_list_to_dotbracket(bp_list, seq_len,offset=0,full_nt_naming=False):
    db = "." * seq_len

    # group into bps that are not intertwined and can use same brackets!
    bp_list = bp_list.copy()
    # if we have chain name in their,add the length of the sequence to each
    if full_nt_naming:
        # if multiple chain we need more
        # this assume homo multimer
        chains = get_chains(bp_list,chains=set())
        db = "." * (seq_len * len(chains))
        bp_list,translation_dict = get_residue_numbers(bp_list,chain_add=seq_len,chain_factors={},res_dict={})

    bp_list = [[nuc+offset for nuc in pair] for pair in bp_list ]
    groups = _group_into_non_conflicting_bp(bp_list)

    # all bp that are not intertwined get (), but all others are
    # groups to be nonconflicting and then asigned (), [], {}, <> by group
    chars_set = [("(", ")"), ("(", ")"), ("[", "]"), ("{", "}"), ("<", ">")]
    alphabet = [(chr(lower), chr(upper))
                for upper, lower in zip(list(range(65, 91)), list(range(97, 123)))]
    chars_set.extend(alphabet)

    if len(groups) > len(chars_set):
        print(f"WARNING: PK too complex, not enough brackets to represent it, {len(groups)}.")
    for group, chars in zip(groups, chars_set):
        for bp in group:
            if bp[0]>bp[1]:
                bpA = bp[1]
                bpB = bp[0]
            else:
                bpB = bp[1]
                bpA = bp[0]
            db = db[:bpA] + chars[0] + \
                db[bpA + 1:bpB] + chars[1] + db[bpB + 1:]
    if full_nt_naming:
        if len(chains)>1:
            # put chain breaks
            db = "*".join([db[i*seq_len:(i+1)*seq_len] for i in range(len(chains))])
    return db
    

def _get_list_bp_conflicts(bp_list):
    '''given a bp_list gives the list of conflicts bp-s which indicate PK structure
    Args:
            bp_list: of list of base pairs where the base pairs are list of indeces of the bp in increasing order (bp[0]<bp[1])
    returns:
            List of conflicting basepairs, where conflicting is pairs of base pairs that are intertwined.
    '''
    if len(bp_list) <= 1:
        return []
    else:
        current_bp = bp_list[0]
        conflicts = []
        for bp in bp_list[1:]:
            if (bp[0] < current_bp[1] and current_bp[1] < bp[1]):
                conflicts.append([current_bp, bp])
        return conflicts + _get_list_bp_conflicts(bp_list[1:])

def _get_non_redudant_bp_list(conflict_list):
    ''' given a conflict list get the list of nonredundant basepairs this list has

    Args:
            conflict_list: list of pairs of base_pairs that are intertwined basepairs
    returns:
            list of basepairs in conflict list without repeats
    '''
    non_redudant_bp_list = []
    for conflict in conflict_list:
        if conflict[0] not in non_redudant_bp_list:
            non_redudant_bp_list.append(conflict[0])
        if conflict[1] not in non_redudant_bp_list:
            non_redudant_bp_list.append(conflict[1])
    return non_redudant_bp_list

def _group_into_non_conflicting_bp(bp_list):
    ''' given a bp_list, group basepairs into groups that do not conflict

    Args
            bp_list: list of base_pairs

    Returns:
            groups of baspairs that are not intertwined
    '''
    conflict_list = _get_list_bp_conflicts(bp_list)

    non_redudant_bp_list = _get_non_redudant_bp_list(conflict_list)
    bp_with_no_conflict = [
        bp for bp in bp_list if bp not in non_redudant_bp_list]
    groups = [bp_with_no_conflict]
    while non_redudant_bp_list != []:
        current_bp = non_redudant_bp_list[0]
        current_bp_conflicts = []
        for conflict in conflict_list:
            if current_bp == conflict[0]:
                current_bp_conflicts.append(conflict[1])
            elif current_bp == conflict[1]:
                current_bp_conflicts.append(conflict[0])
        max_group = [
            bp for bp in non_redudant_bp_list if bp not in current_bp_conflicts]
        to_remove = []
        for i, bpA in enumerate(max_group):
            for bpB in max_group[i:]:
                if bpA not in to_remove and bpB not in to_remove:
                    if [bpA, bpB] in conflict_list or [bpB, bpA] in conflict_list:
                        to_remove.append(bpB)
        group = [bp for bp in max_group if bp not in to_remove]
        groups.append(group)
        # remove group from list
        non_redudant_bp_list = [bp for bp in non_redudant_bp_list if bp not in group] # current_bp_conflicts
        # conflict_list = [conflict for conflict in conflict_list if conflict[0] not in group and conflict[1] not in group]
    return groups

def get_bp_lists(struct_list,target,single_struct=False,offset=0):
    bp_lists = []
    singlets = []
    cross_pair_lists = []
    if single_struct:
        struct_list = [struct_list]
        
    for struct in struct_list:
        singlet_list = []
        bp_list = convert_dotbracket_to_bp_list(struct, allow_pseudoknots=True,silent_values=['-','&'])

        # remove any pair with somethin unresolved
        bp_list = [pair for pair in bp_list if is_resolved(pair[0],target) and is_resolved(pair[1],target)]
        # not any(num in pair for num in ignore_list)]

        bp_list = [[nuc+offset for nuc in pair] for pair in bp_list ]
        helices = get_helices(bp_list)
        for helix in helices:
            if len(helix)==1:
                singlet_list.append(helix[0])

        groups = _group_into_non_conflicting_bp(bp_list)
        # bp_list_no_pk = groups[0]
        bp_list_pk = [bp for group in groups[1:] for bp in group]
        cross_pair_lists.append(bp_list_pk)
        bp_lists.append(bp_list)
        singlets.append(singlet_list)
    if single_struct:
        return bp_lists[0],singlets[0],cross_pair_lists[0]
    else:
        return bp_lists,singlets,cross_pair_lists

def get_crossed_singlets(bp_lists,target,single_struct=False,nt_num_only=True,verbose=False):
    singlets = []
    cross_pair_lists = []

    if single_struct:
        bp_lists = [bp_lists]
        
    for bp_list in bp_lists:
        singlet_list = []
        
        bp_list_res,translation_dict = get_residue_numbers(bp_list,chain_add=10000,chain_factors={},res_dict={})

        #problem_bp_list = [pair for pair in bp_list_res if any(num in pair for num in ignore_list)]
        problem_bp_list = [pair for pair in bp_list_res if not (is_resolved(pair[0]%10000,target) and is_resolved(pair[1]%10000,target))]
        if verbose:
            if len(problem_bp_list)>0: print("unresolved in bp_list", problem_bp_list)
        helices = get_helices(bp_list,nt_num_only=nt_num_only)
        for helix in helices:
            if len(helix)==1:
                singlet_list.append(helix[0])
                
        groups = _group_into_non_conflicting_bp(bp_list_res)
        
        bp_list_pk = [[translation_dict[nt] for nt in bp] for group in groups[1:] for bp in group]
        
        cross_pair_lists.append(bp_list_pk)
        singlets.append(singlet_list)
    if single_struct:
        return singlets[0],cross_pair_lists[0]
    else:
        return singlets,cross_pair_lists

def get_helices(bp_list, allowed_buldge_len=0,nt_num_only=True):
    bp_list = bp_list.copy()
    if not nt_num_only:
        bp_list,res_dict = get_residue_numbers(bp_list,chain_factors={},res_dict={})
    helices = []
    current_helix = []
    current_helix_return = []
    while bp_list != []:
        current_bp = bp_list.pop(0)
        if current_helix == []:
            current_helix.append(current_bp)
            if not nt_num_only:
                current_helix_return = [[res_dict[current_bp[0]], res_dict[current_bp[1]]]]
        else:
            in_helix_left = list(range(current_helix[-1][0] + 1, current_helix[-1][0] + allowed_buldge_len + 2))
            in_helix_right = list(range(current_helix[-1][1] - allowed_buldge_len - 1, current_helix[-1][1]))
            if current_bp[0] in in_helix_left and current_bp[1] in in_helix_right:
                current_helix.append(current_bp)
                if not nt_num_only:
                    current_helix_return.append([res_dict[current_bp[0]], res_dict[current_bp[1]]])
            else:
                if nt_num_only:
                    helices.append(current_helix)
                else:
                    helices.append(current_helix_return)
                current_helix = [current_bp]
                if not nt_num_only:
                    current_helix_return = [[res_dict[current_bp[0]], res_dict[current_bp[1]]]]
    if nt_num_only:
        helices.append(current_helix)
    else:
        helices.append(current_helix_return)
    return helices


###############################################################################
# F1 scoring for motifs
def remove_uresolved(preds,target,chain_unresolved_same=True):
    # ignore_list

    # for chains we assume we remove the same ignore_list from all
    # obviously only works with the homo-multimers

    #ignore_set = set(ignore_list)
    if chain_unresolved_same:
        # we assume nubmering is correct, ignore residue_name
        ok_res = []
        residues,res_dict = get_residue_numbers(preds,chain_factors={},res_dict={})
        for bp in residues:            
            if all(is_resolved(element, target) for element in bp):
            #is_resolved(bp[0],target) and is_resolved(bp[1],target):
            #if not set(bp).intersection(ignore_set):
                ok_res.append([res_dict[elem] for elem in bp])
        return ok_res
    else:
        print('ERROR: remove_uresolved not yet chain aware implementation, assumes all chains the same')
        return None
    
    
def compare_set_lists(refs,pred,target,zero_division = 0, # ignore_list=[]
                     both_none_value = np.nan, chain_unresolved_same=True,
                    return_recall=False,return_precision=False ):
    
    if not chain_unresolved_same:
        print('ERROR: compare_set_lists not yet chain aware implementation, assumes all chains the same')
        return None
    # MCC = (TP * TN — FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    # ignore list has residues that were unresovled in native
    # we ignore any predictions invovling these
    pred = remove_uresolved(pred,target,chain_unresolved_same=True)

    pred = [get_residue_numbers(x,chain_factors={},res_dict={})[0] for x in pred]
    
    
    f1s = []
    # go through each native structure
    for ref in refs:
        # if both do not have the motif type
        if len(pred) == 0 and len(ref) == 0:
            f1s.append(both_none_value)
        else:
            # assumes things are ordered similiarly
            ref = [get_residue_numbers(x,chain_factors={},res_dict={})[0] for x in ref]
            ref = set(tuple(l) for l in ref)
            pred = set(tuple(l) for l in pred)
        
            TP = len(ref.intersection(pred))
            P = len(ref)
            FP = len(pred - ref)
            FN = len(ref - pred)
        
            recall = TP/P if (P) >0 else zero_division
            precision = TP/(TP+FP) if (TP+FP) >0 else zero_division
            F1 = 2* ((precision*recall)/(precision+recall))  if (precision+recall) > 0 else zero_division
            
            if return_recall:
                f1s.append(recall)
            elif return_precision:
                f1s.append(precision)
            else:
                f1s.append(F1)
    return f1s



###############################################################################
# processing residue nomenclature
def get_residue_numbers(residue_list,chain_add=None,chain_factors={},res_dict={}):
    if len(chain_factors.keys()) == 0 and chain_add is not None:
        # we jsut started so get all chain adders
        # get all chains
        chains = list(get_chains(residue_list,chains=set()))
        for i,chain in enumerate(chains):
            chain_factors[chain] = i*chain_add

    if isinstance(residue_list, list):
        residue_list = residue_list.copy()
        new_list = []
        for item in residue_list:
            new_res, res_dict = get_residue_numbers(item,chain_add=chain_add,chain_factors=chain_factors,res_dict=res_dict)
            new_list.append(new_res)
        return new_list,res_dict

    elif isinstance(residue_list,int):
        # this is just for when we put a no chain information list in this function
        res_dict[residue_list] = residue_list
        return residue_list,res_dict
    else:
        new_res = int(residue_list.split(':')[-1][1:]) 
        if chain_add is not None:
            new_res += chain_factors[get_chain_from_nt(residue_list)]
        res_dict[new_res] = residue_list
        return new_res, res_dict


def get_chains(residue_list, chains):
    if isinstance(residue_list, list):
        for res in residue_list:
            chains = get_chains(res, chains) 
    else:
        chain_identifier = get_chain_from_nt(residue_list)  # Take everything before the last colon
        chains.add(chain_identifier) 
    return chains

def get_residue_from_rosetta(res):
    res_name,chain_res_num = res.split('-')
    chain_res_num = chain_res_num.split(':')
    res_num = chain_res_num[-1]
    chain = ':'.join(chain_res_num[:-1])
    return f'{chain}:{res_name}{res_num}'
    
def get_residues_from_rosetta(residues):
    #print(residues)
    return [get_residue_from_rosetta(res) for res in residues.split('_')]

def get_chain_from_nt(nt):
    return ':'.join(nt.split(':')[:-1])


###############################################################################
# casp ranking
###############################################################################



def reduce_df(df, score_to_choice_best=None,
              static_columns=["target", "gr_code", "model"],
              metric_dict=METRICS,
              participipation_cutoff=None,participation_score=None,
              participation_targets=None,participation_fillna=False):
    ''' 
    Reduce df to "best" score per X 
        (vs multiple scores per X if there are multiple Y)
    Eg for conformations:
    If score_to_choice_best is None then, for each score, the best value over all conformations is selected
    If it is specified then first the best conformation according to that score is selected,
    then all scores are taken from the conformation.

    Args:
        df (DataFrame): DataFrame with all columns in stat_columns and metric_dict, will combine all rows for any other columns
        score_to_choice_best (str): if None then for each score in metric_dict take best value according to agg function, 
            if it is specified then take all values from the same row which is the best score according to this column (default None)
        static_columns (list of str): The columns to group by
        metric_dict (dict str(metric):str(agg function)): the dataframe column names and the function to aggregate scores (eg 'clashscore':'min')

    Returns:
        DataFrame reduced to one value per static_columns by the aggregation in metric_dict
    '''
    columns_combining = [x for x in df.columns if x not in list(
        metric_dict.keys()) + static_columns]
    print("combing the following columns:", columns_combining)

    if score_to_choice_best is not None:
        if metric_dict[score_to_choice_best] == "max":
            best_index = df.groupby(static_columns)[
                score_to_choice_best].idxmax().to_numpy(copy=True)
        elif metric_dict[score_to_choice_best] == "min":
            best_index = df.groupby(static_columns)[
                score_to_choice_best].idxmin().to_numpy(copy=True)
        reduced_df = df.iloc[best_index].drop(columns_combining, axis=1)
    else:
        reduced_df = df.groupby(static_columns).agg(metric_dict).reset_index()

    if participipation_cutoff:

        df_aggregated= reduced_df[reduced_df.target.isin(participation_targets)].copy()
        df_aggregated = df_aggregated.groupby(static_columns[:2]).agg(metric_dict).reset_index()
        if participation_fillna:
            df_aggregated[participation_score] = df_aggregated[participation_score].fillna(1)
        pivot_df = df_aggregated.pivot(index=static_columns[1], 
            columns=static_columns[0], values=participation_score)
        mask = pivot_df.isnull()

        # TARGETS_CHOOSE_BEST is a list of list of targets
        # for every list of targets, if the 2 or more of the target are in the table
        # then I should take whether any of them are not null 
        # then collpase the columns so there is only one column representing each list of targets
        
        new_participation_columns, columns_to_drop= [],[]
        for target_group in TARGETS_CHOOSE_BEST:  
            # Filter the target group to include only targets that are actually columns in pivot_df
            valid_targets = [target for target in target_group if target in pivot_df.columns]

            # If no valid targets are found, skip this group.  This avoids errors.
            if not valid_targets:
                continue # Go to the next target_group

            group_mask = pivot_df.loc[:, valid_targets].isnull()

            # Aggregate participation within the group by checking if ANY valid target is non-null
            participation_within_group = group_mask.all(axis=1)  # Use .any() across columns (targets)
            new_participation_columns.append(participation_within_group.rename(f"participation_{'_'.join(target_group)}"))
            columns_to_drop.extend(valid_targets)
        # Merge new participation columns into pivot_df
        if new_participation_columns:
            merged_participation = pd.concat(new_participation_columns, axis=1)
            mask = pd.concat([mask, merged_participation], axis=1)
            mask = mask.drop(columns=columns_to_drop, errors='ignore')

        # Calculate the overall participation, including the new columns if available.
        # This part now uses the modified or original pivot_df
        # mask is now calculated on the *entire* pivot_df
        #mask = pivot_df.isnull()
        participation = 1 - (mask.sum(axis=1) / len(mask.columns))
        print(participation)
        groups_to_drop = pivot_df[participation <= PARTICIPATION_RATE].index
        reduced_df = reduced_df[~reduced_df[static_columns[1]].isin(groups_to_drop)]
        print('dropping',groups_to_drop)

    return reduced_df




def calc_z(arr, raw):
    '''
    Calculate Z of second array with mean and std of the first

    Args:
        arr (array): is values to calculate mean and std of
        raw (array): is values to calculate z with calculated mean and std

    Returns:
        array of zscore
    '''
    u = np.nanmean(arr)
    s = np.nanstd(arr)
    z = (raw - u) / s
    return z



def get_zscore(arr, negative=False, tolerance_threshold=-2, 
               penalty_threshold=None, mask=None,
              num_filtering_rounds=1):
    '''
    Calculate Zscore of an array

    1. Calculate zscores from the raw scores for all models
    2. Remove outliers - models with zscores below the tolerance threshold (set to -2.0);
    3. Recalculate zscores on the reduced dataset
    4. Assign z-scores below the penalty threshold to the value of this threshold.

    Args:
        arr (array): the array of scores to get z-score of, will ignore nans
        negative (bool): true if the minimum score is the best score (default False)
        threshold (float): after initial mean and std calculations, any Zscore below -2 is ignore for the final mean and
            threshold calculation, the final Zscores are also capped at threshold (default -2)
        mask (array): will ignore any scores in arr where mask is False (default None)

    Returns:
        Array of the Z-scores
    '''

    # copy raw scores to calculate final z later
    raw = arr.astype(float)
    old_raw = np.copy(raw)

    # mask out values for mean and std calcualtion
    if mask is not None:
        raw[~mask] = np.nan

    # set multiple if large or small is good
    if negative:
        mult = -1.
    else:
        mult = 1.

    # calculate z score
    zscore = calc_z(raw * mult, raw * mult)
    # ignore negative outliers and recalculate
    # mean and std, calculate z over all raw values
    for i in range(num_filtering_rounds):
        raw[zscore < tolerance_threshold] = np.nan
        zscore = calc_z(raw * mult, old_raw * mult)
    # cap z score at threshold
    if penalty_threshold is not None:
        zscore[zscore < penalty_threshold] = penalty_threshold
    return zscore


def get_weighted_sum_z(data, weights, prefix=''):
    '''
    Obtain the weight sum of scores
    Weights should sum to 1
    Note does to explicitly require weights to sum to 1

    Args:
        data (DataFrame): contains columns for each prefix+score in weights
        weights (dict str(score):float(weight)): the weights for each score
        prefix (str): if column names have a prefix (eg 'Z_') infront of the 
            scores specified in weights, specify this here (default '')

    Returns:
        Series of the summed scores
    '''
    final_z = None
    for score, weight in weights.items():
        if final_z is None:
            final_z = weight * data[prefix + score].to_numpy()
        else:
            final_z += weight * data[prefix + score].to_numpy()
    return final_z



def get_group_score(df_x, agg="sum>0", score="Z_rna",targets_to_choose_best=TARGETS_CHOOSE_BEST,
    group_col='gr_code',fill_nan=np.nan,fill_not_present=None,target_col='target',weights=[]):
    '''
    Aggregate the scores of each group in gr_code

    Args:
        df (DataFrame): contain gr_code and score (argument specified) column
        agg (str): how to aggregate scores options: 'sum>0','mean','sum' (default 'sum>0')
        score (str): what column to aggregate

    Returns:
        DataFrame of groups and their summed score
    '''
    df = df_x.copy()

    original_nan_pairs = set(df[df[score].isnull()][[group_col, target_col]].copy().itertuples(index=False))

    # choose best score if needed
    for targ_list in targets_to_choose_best:
        avail_targ = [t for t in targ_list if t in df[target_col].unique()]
        if len(avail_targ)<2:
            continue
        best_scores = df[df[target_col].isin(targ_list)].groupby(group_col)[score].apply(_max_with_nan_per_column).to_dict()
        df[score] = df.apply(lambda x: x[score] if x[target_col]!=targ_list[0] else best_scores.get(x[group_col],0), axis=1)
        # remove other rows
        df = df[~df[target_col].isin(targ_list[1:])]

    if fill_not_present is not None:
        # Get all unique group and target combinations
        all_groups = df[group_col].unique()
        all_targets = df[target_col].unique()

        # Create a complete grid of group and target combinations
        complete_grid = pd.DataFrame([(g, t) for g in all_groups for t in all_targets],
                                     columns=[group_col, target_col])

        # Merge with the original DataFrame to find missing combinations
        merged_df = pd.merge(complete_grid, df, on=[group_col, target_col], how='left')

        # Fill missing scores
        merged_df[score] = merged_df[score].fillna(fill_not_present)

        df = merged_df 
    df[score] = df.apply(lambda row: fill_nan if (row[group_col], row[target_col]) in original_nan_pairs else row[score], axis=1)
    vals = df.groupby(group_col)[score]
    if agg == "sum>0":
        return vals.apply(lambda col: col[col > 0].sum()).reset_index()
    elif agg == "mean":
        return vals.mean().reset_index()
    elif agg == "sum":
        return vals.sum().reset_index()

def get_group_score_bootstrap(df_x, agg="sum>0", score="Z_rna", targets_to_choose_best=TARGETS_CHOOSE_BEST,
    group_col='gr_code', fill_nan=np.nan, fill_not_present=None, n_bootstraps=100,target_col='target',weights=[]):

    '''
    Aggregate the scores of each group in gr_code, with bootstrapping.

    Args:
        df_x (DataFrame): DataFrame containing gr_code, target, and score columns.
        agg (str): How to aggregate scores. Options: 'sum>0', 'mean', 'sum'. (default 'sum>0')
        score (str): The name of the score column.
        targets_to_choose_best (list of lists):  Targets to choose the best score for within each inner list.
        group_col (str): The name of the group column. (default 'gr_code')
        fill_nan (numeric): Value to fill NaN scores with. (default np.nan)
        fill_not_present (numeric): Value to fill missing group/target combinations with. (default None)
        n_bootstraps (int): Number of bootstrap samples. (default 100)

    Returns:
        DataFrame:  A DataFrame containing the group scores, and a bootstrap_variances dictionary. The dictionary contains the bootstrapped variance for each group.
    '''
    df = df_x.copy()

    # Get unique targets for bootstrapping
    targets = df[target_col].unique()

    # Store the results from the original function
    original_results = get_group_score(df_x, agg=agg, score=score, targets_to_choose_best=targets_to_choose_best,
                            group_col=group_col, fill_nan=fill_nan, fill_not_present=fill_not_present,target_col=target_col,weights=weights)
    original_results = original_results.set_index(group_col)

    # Bootstrap loop
    #for group in original_results.index:
    bootstrap_samples = pd.DataFrame(index=original_results.index)#[]
    bootstrap_data= {}

    for n in tqdm(range(n_bootstraps)):
        # Bootstrap targets (sample with replacement)
        bootstrapped_targets = np.random.choice(targets, size=len(targets), replace=True)
        df_bootstrapped = []
        for t in bootstrapped_targets:
            df_bootstrapped.append(df_x[df_x[target_col]==t])
        df_bootstrapped = pd.concat(df_bootstrapped)
        #print('input',df_bootstrapped)
        # Recalculate scores using the bootstrapped data
        bootstrapped_scores = get_group_score(df_bootstrapped, agg=agg, score=score, targets_to_choose_best=targets_to_choose_best,
                              group_col=group_col, fill_nan=fill_nan, fill_not_present=fill_not_present,target_col=target_col)

        # Store the result from this bootstrap run
        bootstrapped_scores = bootstrapped_scores.set_index(group_col)
        bootstrap_data[n] = bootstrapped_scores[score] 
    bootstrap_samples = pd.DataFrame.from_dict(bootstrap_data, orient='columns') # Now make the dataframe from the dict


    # Calculate the variance for the group across the bootstrap samples
    bootstrap_samples[f'mean_{score}'] = bootstrap_samples[list(range(n_bootstraps))].mean(axis=1)
    bootstrap_samples[f'sem_{score}'] = bootstrap_samples[list(range(n_bootstraps))].sem(axis=1)
    bootstrap_samples[f'std_{score}'] = bootstrap_samples[list(range(n_bootstraps))].std(axis=1)
    bootstrap_samples[f'q2.5_{score}'] = bootstrap_samples[list(range(n_bootstraps))].quantile(0.025,axis=1) 
    bootstrap_samples[f'q97.5_{score}'] = bootstrap_samples[list(range(n_bootstraps))].quantile(0.975,axis=1)  
    bootstrap_samples[f'q20_{score}'] = bootstrap_samples[list(range(n_bootstraps))].quantile(0.2,axis=1) 
    bootstrap_samples[f'q80_{score}'] = bootstrap_samples[list(range(n_bootstraps))].quantile(0.8,axis=1)  
    bootstrap_samples[f'q15.9_{score}'] = bootstrap_samples[list(range(n_bootstraps))].quantile(0.159,axis=1) 
    bootstrap_samples[f'q84.1_{score}'] = bootstrap_samples[list(range(n_bootstraps))].quantile(0.841,axis=1)  

    return original_results, bootstrap_samples


###############################################################################
# plotting code

# code to make the GDT over all targets and graph plot
def overview_performance_clustermap(plot_df,x_axis,y_axis,value,mid_value,
    participation_rate,score_name,target_colors,group_colors,save,rotate_x=90,label_size=8,
    figsize=(8, 7),scores=None):
    # note all these should be max metrics
    df_aggregated = plot_df.groupby([y_axis, x_axis])[value].max().reset_index()
    pivot_df = df_aggregated.pivot(index=y_axis, columns=x_axis, values=value)
    
    # make a color mpa
    cmap = mcolors.LinearSegmentedColormap.from_list("", ["red", "white", "blue"])
    cmap.set_bad('gray', alpha=1.0) #set NaN to grey
    mask = pivot_df.isnull()
    
    # clustermap algorithm of clsutering cannot deal with nans so set value
    # and reduce to only those participating in most
    #pivot_df[mask] = 0
    participation = 1-(mask.sum(axis=1) / len(mask.columns))
    groups_to_drop = pivot_df[participation <= participation_rate].index
    pivot_df = pivot_df.drop(groups_to_drop)
    mask = mask.drop(groups_to_drop)
    pivot_df.fillna(0, inplace=True)

    # remove R0 targets
    columns_to_keep = [col for col in pivot_df.columns if not col.startswith('R0')]
    pivot_df = pivot_df[columns_to_keep]
    mask = mask[columns_to_keep]

    # If a target (pivot_df.columns) is not in target_colors add a row for that target with all "grey"
    missing_targets = set(pivot_df.columns) - set(target_colors.index)
    for target in missing_targets:
        new_row = pd.DataFrame({col: ['lightgray'] for col in target_colors.columns}, index=[target])
        target_colors = pd.concat([target_colors, new_row])

    # Calculate the row sums and reorder the pivot table
    row_sums = pivot_df.sum(axis=1)
    pivot_df = pivot_df.loc[row_sums.sort_values(ascending=False).index]
    mask = mask.loc[row_sums.sort_values(ascending=False).index]

    g = sns.clustermap(pivot_df, 
                       cmap=cmap, 
                       center=mid_value, 
                        vmin=0,vmax=1, # make cmap is all the way from 0 to 1
                       mask=mask, # grey out now predicted
                       metric="euclidean", 
                       figsize=figsize,
                      #colors_ratio=(0.8,0.8), # make it square
                      yticklabels=True,
                      xticklabels=True, # make sure all row/column labels shows
                       col_colors = target_colors,
                       row_colors = group_colors,
                       colors_ratio=0.02,
                      )
    if scores is not None:        
        linkage_matrix = g.dendrogram_row.linkage

        def order_leaves_by_scores(linkage_matrix, scores, N=3):
            from scipy.cluster.hierarchy import to_tree
            tree, node_dict = to_tree(linkage_matrix, rd=True)
        
            def get_order(node):
                """
                Function to get the leaf order starting from the highest score.
                """
                if node.is_leaf():
                    return [node.id], [scores[node.id]]
                
                left_order, left_scores = get_order(node.get_left())
                right_order, right_scores = get_order(node.get_right())
                
                left_scores_sorted = sorted(left_scores, reverse=True)[:N]
                right_scores_sorted = sorted(right_scores, reverse=True)[:N]
                
                left_highest = sum(left_scores_sorted)/len(left_scores_sorted)
                right_highest = sum(right_scores_sorted)/len(right_scores_sorted)
                
                if left_highest >= right_highest:
                    return left_order + right_order, left_scores + right_scores
                else:
                    return right_order + left_order, right_scores + left_scores
        
            ordered_ids, _ = get_order(tree)
            return ordered_ids

        leaf_order = order_leaves_by_scores(linkage_matrix, scores.loc[pivot_df.index], N=5)

        plt.clf()
        g = sns.clustermap(pivot_df.iloc[leaf_order], 
                           cmap=cmap, 
                           center=mid_value, 
                           vmin=0, vmax=1,  # Make cmap full scale
                           mask=mask.iloc[leaf_order],  # Grey out not predicted
                           metric="euclidean", method='average',
                           figsize=figsize,
                           yticklabels=True,
                           xticklabels=True,  # Ensure all labels show
                           col_colors=target_colors,
                           row_colors=group_colors.loc[pivot_df.iloc[leaf_order].index],
                           #colors_ratio=0.02,#row_linkage=new_linkage_matrix,
                           row_cluster=False,
                           colors_ratio=(0.01, 0.02),
                          dendrogram_ratio=0.1)

    g.ax_heatmap.tick_params(axis='both', labelsize=label_size) 
    g.ax_heatmap.set_ylabel("Group", fontsize=12) 
    g.ax_heatmap.set_xlabel("Target", fontsize=12) 
    for tick in g.ax_heatmap.get_xticklabels():
        tick.set_rotation(rotate_x)
        tick.set_ha("right")
    
    # colormap edits
    cax = plt.gcf().axes[-1] # Get the last axes (colorbar)
    tick_locator = ticker.FixedLocator([0, mid_value, 1])
    cax.yaxis.set_major_locator(tick_locator)
    pos = cax.get_position()
    cax.set_position([pos.x0 + 0.05, pos.y0+0.04, pos.width * 0.6, pos.height*0.6]) 
    cax.tick_params(labelsize=8)
    plt.title(score_name, loc='left', fontsize=12)
    plt.savefig(f"{save}.png",dpi=400, bbox_inches='tight', transparent=True)
    plt.savefig(f"{save}.svg",dpi=400, bbox_inches='tight', transparent=True)
    return leaf_order, linkage_matrix, pivot_df.iloc[leaf_order]

def plot_ss_accuracy(preds_top1,score):

    fig = plt.figure(figsize=(10,8))
    gs = GridSpec(3, 5,wspace=1)
    ax_main = plt.subplot(gs[0:3, :4])
    ax_yDist = plt.subplot(gs[0:3, 4:], sharey=ax_main)
    
    rank_order = preds_top1.groupby('group').sum().reset_index().sort_values(score, ascending=False).group.values
    preds_top1['gr_code'] = pd.Categorical(preds_top1['group'],ordered=True,categories=rank_order)
    preds_pivot = preds_top1.pivot_table(index="gr_code",columns="target",values=score, dropna=False) # f1_singlet
    mask = np.isnan(preds_pivot)
    g = sns.heatmap(preds_pivot,cmap="seismic_r",yticklabels=True,xticklabels=True,
                   cbar_kws={'shrink': 0.6},mask=mask,ax=ax_main)
    t=g.set_yticklabels(labels=g.get_yticklabels(),rotation = 0,size=6)
    t=g.set_xticklabels(labels=g.get_xticklabels(),rotation = 90,size=6)
    g.set_facecolor('lightgrey')

    plt.barh(g.get_yticks(),width=preds_top1.groupby('group').sum().reset_index().sort_values(score, ascending=False)[score],
             color="grey")
    ax_yDist.set_xlabel("Sum")
    ax_yDist.set_xlim(0,len(preds_pivot.columns))
    t=ax_yDist.set_yticklabels(labels=g.get_yticklabels(),rotation = 0,size=6)
    #ax_yDist.spines['right'].set_visible(False)
    ax_yDist.spines['top'].set_visible(False)


def plot_barplot_SS(preds_top1, scores, palette,figsize=(10, 6),top_N=None,y_highlight=[],save=None,
                ticklabel_size=5,title='',yrot=0,rank_type='sum',targets_participating=None,
                group_col='group',targets_to_choose_best=TARGETS_CHOOSE_BEST):
    df_copy = preds_top1.copy()

    if rank_type == 'sum':
        group_score = []
        for score in scores:
            group_score.append(get_group_score(df_copy, agg="sum", score=score,
                    targets_to_choose_best=targets_to_choose_best,group_col=group_col,
                   fill_nan=1,fill_not_present=0).sort_values(score,ascending=False).set_index(group_col))
        group_score = pd.concat(group_score,axis=1)

    elif rank_type == 'mean':
        group_score = []
        for score in scores:
            group_score.append(get_group_score(df_copy, agg="mean", score=score, # mean ignores nans so if 6 values, 2 nan, will noramlize sum by 4
                    targets_to_choose_best=targets_to_choose_best,group_col=group_col,
                   fill_nan=np.nan,fill_not_present=0).sort_values(score,ascending=False).set_index(group_col))
        group_score = pd.concat(group_score,axis=1)
        
    # summary ranking just sum, no normalization TODO ???
    rank_order = group_score.sum(1)
    rank_order = rank_order.index.values

    if not top_N is None:
        rank_order = rank_order[:top_N]
    preds_long = group_score.reset_index().melt(id_vars=['group'], value_vars=scores, 
        var_name='score_type', value_name='score')
    
    fig, ax = plt.subplots(figsize=figsize)  
    sns.barplot(preds_long,x='score',y='group',order=rank_order,
                hue='score_type',hue_order=scores,palette=palette,legend=False,
                width=0.8,linewidth=0.5, edgecolor="black")
    ax.spines[['top','right']].set_visible(False)
    #ax.set_xlim(0,1)
    ax.set_title(title)
    t=ax.set_yticklabels(labels=ax.get_yticklabels(),rotation = yrot,size=ticklabel_size)
    t=ax.set_xticklabels(labels=ax.get_xticklabels(),rotation = 0,size=ticklabel_size+1)
    ax.set_ylabel("Group")
    ax.set_xlabel("Sum F1-score")
    for label in ax.get_yticklabels():
        if label.get_text() in y_highlight:
            label.set_color('grey')
    if save:
        plt.savefig(f"{save}.svg",dpi=400, bbox_inches='tight', transparent=True)
        plt.savefig(f"{save}.png",dpi=400, bbox_inches='tight', transparent=True)

def calc_ci_limits(score,min_boundary,max_boundary,nonneg=True):
    lower_bound = score-min_boundary
    upper_bound = max_boundary-score
    if nonneg:
        lower_bound = [max(0,x) for x in lower_bound]
        upper_bound = [max(0,x) for x in upper_bound]
    return [lower_bound,upper_bound]

def plot_heat_map(df,score,figsize=(12,8),cmap='seismic_r',
    savefig=None,targets_to_choose_best=TARGETS_CHOOSE_BEST,TS=False,y_text_size=6,
    h_bar_color='grey',rank_order=None,
    rank_type=None,group_col='gr_code',
    vmax=None,vmin=None,
    y_highlight=[],cbar_label='',ticklabel_size=6,
    num_bootstrap=0,
    split_dims = (52,76),
    target_col='target',target_order=None):


    # setup figure
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(3, 100)
    ax_main = plt.subplot(gs[0:3, :split_dims[0]])
    ax_yDist = plt.subplot(gs[0:3, split_dims[1]:], sharey=ax_main)

    # get ranks
    df_copy = df.copy()
    if "Z_" == score[:2]:
        group_score,score_boot = get_group_score_bootstrap(df_copy,score=score,
            targets_to_choose_best=targets_to_choose_best,group_col=group_col,
            n_bootstraps=num_bootstrap,target_col=target_col)
    elif rank_type == 'sum':
        group_score,score_boot = get_group_score_bootstrap(df_copy, agg="sum", score=score,
                    targets_to_choose_best=targets_to_choose_best,group_col=group_col,
                   fill_nan=1,fill_not_present=0,
                   n_bootstraps=num_bootstrap,target_col=target_col)
    elif rank_type == 'mean':
        group_score,score_boot = get_group_score_bootstrap(df_copy, agg="mean", score=score, # mean ignores nans so if 6 values, 2 nan, will noramlize sum by 4
                    targets_to_choose_best=targets_to_choose_best,group_col=group_col,
                   fill_nan=np.nan,fill_not_present=0,
                   n_bootstraps=num_bootstrap,target_col=target_col)
    else:
        group_score,score_boot = get_group_score_bootstrap(df_copy,score=score,agg='sum',
            targets_to_choose_best=targets_to_choose_best,
                   n_bootstraps=num_bootstrap,target_col=target_col)

    group_score = group_score.reset_index().sort_values(score, ascending=False)
    if rank_order is None:
        rank_order = group_score[group_col].values
        
    # get an ordered pivot table with all scores
    df_copy[group_col] = pd.Categorical(df_copy[group_col],ordered=True,categories=rank_order)
    score_boot= score_boot.reindex(rank_order)

    # no need for fill_value anymore
    nan_marker = -99999  
    df_copy[score] = df_copy[score].fillna(nan_marker)
    pivot_all = df_copy.pivot_table(index=group_col, columns=target_col, values=score,
                                    dropna=False, fill_value=np.nan)

   
    

    # plot the heat map centered on passing score if applicable
    if score in METRIC_PASS_SCORE:
        center = METRIC_PASS_SCORE[score]
    else:
        center = 0
    if target_order is not None:
        pivot_all = pivot_all[target_order]
    
    # create anottation to mark out the ones that are nan before
    nan_mask = (pivot_all == nan_marker) 
    pivot_all = pivot_all.replace(nan_marker, np.nan)
    # mask out any nonparitipating and originally nans
    mask = np.isnan(pivot_all)
    g = sns.heatmap(pivot_all,cmap=cmap,ax=ax_main, center=center,
                   yticklabels=True,xticklabels=True,#[f"{'TS' if TS else ''}{'0'*(3-len(str(x)))}{x}" for x in rank_order], # {group_names[x]}:( TODO when have names
                   cbar_kws={'shrink': 0.6,'location':'top','label':score},
                   mask=mask,vmin=vmin,vmax=vmax,)
                    #annot=nan_locs, fmt='',annot_kws={"fontsize":15}) # annot=True, 
    # Vectorized annotation using np.where
    row_indices, col_indices = np.where(nan_mask)  # Find row, col indices where nan_mask is True
    for row, col in zip(row_indices, col_indices):
        ax_main.text(col + 0.5, row + 0.5, '-', ha='center', va='center', color='dimgrey', fontsize=4)  # Annotate 'o'

    g.set_xticks([x+0.5 for x in range(len(pivot_all.columns))])
    if target_order is not None:
        g.set_xticklabels(target_order, rotation=90, ha='center',size=ticklabel_size)
    else:
        g.set_xticklabels(pivot_all.columns, rotation=90, ha='center',size=ticklabel_size)
    g.set_xlabel('Target')
    t=g.set_yticklabels(labels=g.get_yticklabels(),rotation = 0,size=y_text_size)
    g.set_ylabel("")
    #g.set_title(score)
    g.set_facecolor('lightgrey')
    g.set_ylabel('Group')
    cbar = g.collections[0].colorbar
    cbar.ax.tick_params(labelsize=ticklabel_size)
    cbar.ax.set_position([0.17, 0.75, 0.3, 0.02])
    cbar.ax.xaxis.set_ticks_position('bottom')
    cbar.ax.set_title(cbar_label, fontsize=ticklabel_size*2) 
    cbar.ax.tick_params(axis='x', length=0.8)  # Use 'length' for tick size

    
    # Move tick labels closer to the colorbar
    for label in cbar.ax.get_xticklabels():
        label.set_y(0.2)  # Adjust the vertical position of the labels closer to the colorbar

    # plot horizantal bars
    #print(group_score.set_index(group_col)[score],score_boot[f'q20_{score}'], score_boot[f'q80_{score}'],group_score[score])
    #print(group_score.set_index(group_col)[score]-score_boot[f'q20_{score}'], score_boot[f'q80_{score}']-group_score[score])
    #ci = [(group_score.set_index(group_col)[score]-score_boot[f'q20_{score}']).to_list(), (score_boot[f'q80_{score}']-group_score.set_index(group_col)[score]).to_list()]
    #print(ci)
    ci = calc_ci_limits(group_score.set_index(group_col).reindex(rank_order)[score],score_boot[f'q15.9_{score}'],score_boot[f'q84.1_{score}'])
    plt.barh(g.get_yticks(),width=group_score.set_index(group_col).reindex(rank_order)[score],color=h_bar_color,xerr=ci) #score_boot[f'sem_{score}'].to_list()
    if "Z_" == score[:2]:
        ax_yDist.set_xlabel("Sum Z>0")
    else:
        ax_yDist.set_xlabel("Sum")

    t=ax_yDist.set_yticklabels(labels=g.get_yticklabels(),rotation = 0,size=y_text_size)
    ax_yDist.spines['right'].set_visible(False)
    ax_yDist.spines['top'].set_visible(False)

    group_score = group_score.set_index(group_col)
    if 'AF3-server - S 304' in group_score.index: #pivot_all.index:
        # Plot AF3 baseline
        x_coord = group_score[score]['AF3-server - S 304'] 
        ax_yDist.axvline(x=x_coord, color='gold', linestyle='--', 
            linewidth=1, label='AF3-server - S 304')


    ax_main_pos = ax_main.get_position()  # Get position of the main axis
    ax_yDist.set_position([ax_yDist.get_position().x0, ax_main_pos.y0, ax_yDist.get_position().width, ax_main_pos.height])

    for ax in [ax_main,ax_yDist]:
        for label in ax.get_yticklabels():
            if label.get_text() in y_highlight:
                label.set_color('grey')
    if savefig is not None:
        plt.savefig(f'{savefig}.png',dpi=400, bbox_inches='tight')
        plt.savefig(f'{savefig}.svg',dpi=400, bbox_inches='tight')

    return rank_order,group_score,ci


def plot_top_n(df,score,targets_to_choose_best,N,ax,tick_label_size=8,
    num_bootstrap=0,group_col='gr_code',h_bar_color='grey',color_by=None):
    df_copy = df.copy()
    
    #group_score = get_group_score(df_copy,score=score,targets_to_choose_best=targets_to_choose_best)
    #group_score["tie"] = get_group_score(df_copy,agg="mean",score=score,targets_to_choose_best=targets_to_choose_best)[score]
    #top_n = group_score.sort_values(score, ascending=False).head(N)
    #rank_order = top_n.gr_code.values

    if "Z_" == score[:2]:
        group_score,score_boot = get_group_score_bootstrap(df_copy,score=score,
            targets_to_choose_best=targets_to_choose_best,group_col=group_col,
            n_bootstraps=num_bootstrap)
    elif rank_type == 'sum':
        group_score,score_boot = get_group_score_bootstrap(df_copy, agg="sum", score=score,
                    targets_to_choose_best=targets_to_choose_best,group_col=group_col,
                   fill_nan=1,fill_not_present=0,
                   n_bootstraps=num_bootstrap)
    elif rank_type == 'mean':
        group_score,score_boot = get_group_score_bootstrap(df_copy, agg="mean", score=score, # mean ignores nans so if 6 values, 2 nan, will noramlize sum by 4
                    targets_to_choose_best=targets_to_choose_best,group_col=group_col,
                   fill_nan=np.nan,fill_not_present=0,
                   n_bootstraps=num_bootstrap)
    else:
        group_score,score_boot = get_group_score_bootstrap(df_copy,score=score,agg='sum',
            targets_to_choose_best=targets_to_choose_best,
                   n_bootstraps=num_bootstrap)
    
    group_score = group_score.reset_index().sort_values(score, ascending=False)
    rank_order = group_score[group_col].values
    score_boot= score_boot.reindex(rank_order)

    #rank_order = group_score.sort_values(score, ascending=True).gr_code.values
    plot_fig = df_copy[df_copy[group_col].isin(rank_order)].copy()

    #print(group_score)
    #print(score_boot)
    ci = calc_ci_limits(group_score.set_index(group_col).reindex(rank_order)[score],score_boot[f'q15.9_{score}'],score_boot[f'q84.1_{score}'])
    ax.barh(rank_order[:N],width=group_score.set_index(group_col).reindex(rank_order)[score][:N],
            color=h_bar_color,xerr=[c[:N] for c in ci]) #score_boot[f'sem_{score}'].to_list()
    if color_by:
        # assume scores is 1D
        print(preds_top1.head())
        #preds_top1_sorted = preds_top1.sort_values(color_by, key=lambda x: pd.Categorical(x, categories=order_color, ordered=True))

        df_cumsum = preds_top1.groupby(['group', color_by])[scores[0]].sum().unstack()[order_color].cumsum(axis=1).stack().reset_index()
        df_cumsum.rename(columns={0: scores[0]}, inplace=True)

        print(df_cumsum.head())
        ax = sns.barplot(df_cumsum,x=scores[0],y='group',order=rank_order,dodge=False,
                hue=color_by,hue_order=order_color[::-1],palette=palette,legend=True,
                linewidth=0.5, edgecolor="black",width=0.6,ax=ax) # ,
        handles, labels = ax.get_legend_handles_labels()
        legend = ax.legend(handles[::-1], labels[::-1], title=None, loc='lower right', fontsize=8)
    #plot_fig[group_col] = pd.Categorical(plot_fig[group_col],ordered=True,categories=rank_order)
    #pivot_all = plot_fig.pivot_table(index=group_col,columns="target",values=score, dropna=False)
    #mask = np.isnan(pivot_all)
    
    #pivot_all["sum z>0"] = pivot_all.where(pivot_all > 0).sum(1)
    #ax.barh(pivot_all.index,width=pivot_all["sum z>0"],color=color,height=height)
    ax.set_xlabel("Sum Z>0",size=tick_label_size*1.2)
    ax.invert_yaxis()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    t=ax.set_yticklabels(labels=ax.get_yticklabels(),size=tick_label_size)
    t=ax.set_xticklabels(labels=ax.get_xticklabels(),size=tick_label_size*1.2)
    group_score = group_score.set_index(group_col)
    if 'AF3-server - S 304' in group_score.index:
        # Plot AF3 baseline
        x_coord = group_score[score]['AF3-server - S 304'] 
        ax.axvline(x=x_coord, color='gold', linestyle='--', linewidth=1, label='AF3-server - S 304')

# OBSOLETE
'''
def reduce_score_list_overlapping_target_df(scores, direction='max',nan_best=True):

    reduced_df = scores.copy()
    for targets_ in MONOMER_TARGETS_CHOOSE_BEST:
        targets = [t for t in targets_ if t in scores.target.unique()]

        if len(targets)<2:
            continue

        filtered_scores = reduced_df[reduced_df['target'].isin(targets)].copy()
        filtered_scores = filtered_scores.drop('target',axis=1)
        if direction == 'max':
            best_scores = filtered_scores.groupby('group').agg(max_with_nan_per_column).reset_index()
        elif direction == 'min':
            best_scores = filtered_scores.groupby('group').agg(min_with_nan_per_column).reset_index()
        # remove rows with target in targets
        reduced_df = reduced_df[~reduced_df['target'].isin(targets)]
        # add in rows with the reduced score with target = '_'.join(targets) 
        best_scores['target'] = '_'.join(targets)  
        reduced_df = pd.concat([reduced_df, best_scores], ignore_index=True)
    return reduced_df


def plot_heatmap_SS(preds_top1,score,figsize=(10,6),cmap="seismic_r",
                   y_highlight=[],
                  cbar_label='',save=None,ticklabel_size=6,rank_type='sum'):

    fig = plt.figure(figsize=figsize)
    preds_top1 = preds_top1.copy()
    preds_top1_reduced = reduce_score_list_overlapping_target_df(preds_top1)

    if rank_type=='sum':
        # rank by sum
        rank_order = preds_top1_reduced.groupby('group').sum().reset_index().sort_values(score, ascending=False).group.values
    elif rank_type=='sum_norm':
        # rank by sum, with nan being a score of 1
        rank_order_sum_nan_one = (preds_top1_reduced.groupby('group')
            .apply(lambda x: x[score].fillna(1).sum())  # Replace NaN with 1 and then sum
            .reset_index(name='total_score').sort_values('total_score', ascending=False))
        rank_order = rank_order_sum_nan_one['group'].values
    elif rank_type=='mean':
        # rank by nan mean (aka sum and then normalie by number of non nan values)
        rank_order_nan_mean = (preds_top1_reduced.groupby('group')
            .apply(lambda x: x[score].mean(skipna=True))
            .reset_index(name='mean_score').sort_values('mean_score', ascending=False)) 
        rank_order = rank_order_nan_mean['group'].values
    
    preds_top1['gr_code'] = pd.Categorical(preds_top1['group'],ordered=True,categories=rank_order)
    preds_pivot = preds_top1.pivot_table(index="gr_code",columns="target",values=score, dropna=False) # f1_singlet
    mask = np.isnan(preds_pivot)
    g = sns.heatmap(preds_pivot,cmap=cmap,yticklabels=True,xticklabels=True,
                   cbar_kws={'shrink': 0.6,'orientation':'horizontal','location':'top'},
                    mask=mask,vmin=0,vmax=1)
    cbar = g.collections[0].colorbar
    cbar.ax.tick_params(labelsize=6, length=1)   
    cbar.ax.set_title(cbar_label, fontsize=8) 
    cbar.ax.xaxis.set_ticks_position('bottom') 
    cbar.ax.xaxis.set_label_position('bottom')
    t=g.set_yticklabels(labels=g.get_yticklabels(),rotation = 0,size=ticklabel_size)
    t=g.set_xticklabels(labels=g.get_xticklabels(),rotation = 90,size=ticklabel_size)
    g.set_facecolor('lightgrey')
    g.set_ylabel('Group')
    g.set_xlabel('Target')
    for label in g.get_yticklabels():
        if label.get_text() in y_highlight:
            label.set_color('grey')
    if save:
        plt.savefig(f"{save}.svg",dpi=400, bbox_inches='tight')#, transparent=True)
        plt.savefig(f"{save}.png",dpi=400, bbox_inches='tight')#, transparent=True)



def plot_heat_map(df,score,savefig=None,targets_to_choose_best=[],figsize=(12,8),TS=False,y_text_size=6,
    h_bar_color='grey',rank_order=None):
    # TODO this always needs to be best score df...

    df_copy = df.copy()
    
    if "Z_" == score[:2]:
        group_score = get_group_score(df_copy,score=score,targets_to_choose_best=targets_to_choose_best)
        group_score["tie"] = get_group_score(df_copy,agg="mean",score=score,targets_to_choose_best=targets_to_choose_best)[score]
    else:
        group_score = get_group_score(df_copy,score=score,agg='sum',targets_to_choose_best=targets_to_choose_best)

    if rank_order is None:
        rank_order = group_score.sort_values(score, ascending=False).gr_code.values
        
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(3, 100)
    ax_main = plt.subplot(gs[0:3, :52])
    ax_yDist = plt.subplot(gs[0:3, 76:], sharey=ax_main)

    df_copy['gr_code'] = pd.Categorical(df_copy['gr_code'],ordered=True,categories=rank_order)
    pivot_all = df_copy.pivot_table(index="gr_code",columns="target",values=score, dropna=False)

    mask = np.isnan(pivot_all)
    if score in METRIC_PASS_SCORE:
        center = METRIC_PASS_SCORE[score]
    else:
        center = 0
    g = sns.heatmap(pivot_all,cmap="seismic_r",ax=ax_main, center=center,
                   yticklabels=[f"{'TS' if TS else ''}{'0'*(3-len(str(x)))}{x}" for x in rank_order], # {group_names[x]}:( TODO when have names
                   cbar_kws={'shrink': 0.6,'location':'top','label':score},mask=mask) # annot=True, 
    g.set_xticks([x+0.5 for x in range(len(pivot_all.columns))])
    g.set_xticklabels(pivot_all.columns, rotation=90, ha='center',size=6)
    g.set_xlabel('Target')
    t=g.set_yticklabels(labels=g.get_yticklabels(),rotation = 0,size=y_text_size)
    g.set_ylabel("")
    #g.set_title(score)
    g.set_facecolor('lightgrey')
    cbar = g.collections[0].colorbar
    cbar.ax.tick_params(labelsize=6)
    cbar.ax.set_position([0.17, 0.75, 0.3, 0.02])
    cbar.ax.xaxis.set_ticks_position('bottom')
    cbar.ax.tick_params(axis='x', length=0.8)  # Use 'length' for tick size
    
    # Move tick labels closer to the colorbar
    for label in cbar.ax.get_xticklabels():
        label.set_y(0.2)  # Adjust the vertical position of the labels closer to the colorbar

    if "Z_" == score[:2]:
        df_copy = df.copy()
        for targ_list in targets_to_choose_best:
            targ_index = [targ for targ in targ_list if targ in pivot_all.columns]
            best_scores = pivot_all[targ_index].max(1)
            pivot_all[targ_list[0]] = best_scores
            pivot_all[targ_list[1:]] = 0 
        pivot_all["sum z>0"] = pivot_all.where(pivot_all > 0).sum(1)
        plt.barh(g.get_yticks(),width=pivot_all["sum z>0"],color=h_bar_color)
        ax_yDist.set_xlabel("Sum Z>0")
    else:
        cbar.set_label("sum",size=8,ha='right',y=-0.1, rotation=0, va='bottom')
        for targ_list in targets_to_choose_best:
            best_scores = pivot_all[targ_list].max(1)
            pivot_all[targ_list[0]] = best_scores
            pivot_all[targ_list[1:]] = 0 
        pivot_all["sum"] = pivot_all.sum(1)
        plt.barh(g.get_yticks(),width=pivot_all["sum"],color=h_bar_color)
        ax_yDist.set_xlabel("sum")
    # TODO if want to plot another
    #zrna_scores = get_group_score(df_copy)
    #plt.barh(g.get_yticks(),width=zrna_scores.Z_rna,edgecolor="orange", fill=False)
    t=ax_yDist.set_yticklabels(labels=g.get_yticklabels(),rotation = 0,size=y_text_size)
    ax_yDist.spines['right'].set_visible(False)
    ax_yDist.spines['top'].set_visible(False)
    
    #print(pivot_all["sum z>0"])
    #from matplotlib.lines import Line2D
    #custom_lines = [Line2D([0], [0], color="green", lw=4),
                    #Line2D([0], [0], color="orange", lw=4)]

    #ax_yDist.legend(custom_lines, ['Z EM', 'Z RNA'],loc='lower right')

    if 'AF3-server - 304' in pivot_all.index:
        # Plot AF3 baseline
        x_coord = pivot_all["sum z>0" if "Z_" in score else "sum"]['AF3-server - 304'] 
        ax_yDist.axvline(x=x_coord, color='gold', linestyle='--', linewidth=1, label='AF3-server - 304')

    ax_main_pos = ax_main.get_position()  # Get position of the main axis
    ax_yDist.set_position([ax_yDist.get_position().x0, ax_main_pos.y0, ax_yDist.get_position().width, ax_main_pos.height])

    if savefig is not None:
        plt.savefig(f'{savefig}.png',dpi=400, bbox_inches='tight')
        plt.savefig(f'{savefig}.svg',dpi=400, bbox_inches='tight')

    return rank_order

def plot_heat_map(df,score,figsize=(12,8),cmap='seismic_r',
    savefig=None,targets_to_choose_best=TARGETS_CHOOSE_BEST,TS=False,y_text_size=6,
    h_bar_color='grey',rank_order=None,
    rank_type=None,group_col='gr_code',
    vmax=None,vmin=None,
    y_highlight=[],cbar_label='',ticklabel_size=6,
    fill_value=np.nan):

    # setup figure
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(3, 100)
    ax_main = plt.subplot(gs[0:3, :52])
    ax_yDist = plt.subplot(gs[0:3, 76:], sharey=ax_main)

    # get ranks
    df_copy = df.copy()
    if "Z_" == score[:2]:
        group_score = get_group_score(df_copy,score=score,
            targets_to_choose_best=targets_to_choose_best)
        group_score["tie"] = get_group_score(df_copy,agg="mean",
            score=score,targets_to_choose_best=targets_to_choose_best)[score]
    elif rank_type == 'sum':
        group_score = get_group_score(df_copy, agg="sum", score=score,
                    targets_to_choose_best=targets_to_choose_best,group_col=group_col,
                   fill_nan=1,fill_not_present=0) 
    elif rank_type == 'mean':
        group_score = get_group_score(df_copy, agg="mean", score=score, # mean ignores nans so if 6 values, 2 nan, will noramlize sum by 4
                    targets_to_choose_best=targets_to_choose_best,group_col=group_col,
                   fill_nan=np.nan,fill_not_present=0) 
    else:
        group_score = get_group_score(df_copy,score=score,agg='sum',targets_to_choose_best=targets_to_choose_best)

    if rank_order is None:
        rank_order = group_score.sort_values(score, ascending=False)[group_col].values
        
    # get an ordered pivot table with all scores
    df_copy[group_col] = pd.Categorical(df_copy[group_col],ordered=True,categories=rank_order)
    pivot_all = df_copy.pivot_table(index=group_col,columns="target",values=score, 
        dropna=False,fill_value=fill_value)

    # mask out any nonparitipating
    mask = np.isnan(pivot_all)

    # TODO
    # if fill_value is not nan then I want to cmap, it is 
    # the inputted cmap from vmin to vmax, but then yellow outside of that 
    # up to fill_value. The cbar should then only display from vmin to vmax

    # plot the heat map centered on passing score if applicable
    if score in METRIC_PASS_SCORE:
        center = METRIC_PASS_SCORE[score]
    else:
        center = 0

    g = sns.heatmap(pivot_all,cmap=cmap,ax=ax_main, center=center,
                   yticklabels=True,xticklabels=True,#[f"{'TS' if TS else ''}{'0'*(3-len(str(x)))}{x}" for x in rank_order], # {group_names[x]}:( TODO when have names
                   cbar_kws={'shrink': 0.6,'location':'top','label':score},
                   mask=mask,vmin=vmin,vmax=vmax) # annot=True, 
    g.set_xticks([x+0.5 for x in range(len(pivot_all.columns))])
    g.set_xticklabels(pivot_all.columns, rotation=90, ha='center',size=ticklabel_size)
    g.set_xlabel('Target')
    t=g.set_yticklabels(labels=g.get_yticklabels(),rotation = 0,size=y_text_size)
    g.set_ylabel("")
    #g.set_title(score)
    g.set_facecolor('lightgrey')
    g.set_ylabel('Group')
    cbar = g.collections[0].colorbar
    cbar.ax.tick_params(labelsize=ticklabel_size)
    cbar.ax.set_position([0.17, 0.75, 0.3, 0.02])
    cbar.ax.xaxis.set_ticks_position('bottom')
    cbar.ax.set_title(cbar_label, fontsize=ticklabel_size*2) 
    cbar.ax.tick_params(axis='x', length=0.8)  # Use 'length' for tick size
    for label in g.get_yticklabels():
        if label.get_text() in y_highlight:
            label.set_color('grey')
    
    # Move tick labels closer to the colorbar
    for label in cbar.ax.get_xticklabels():
        label.set_y(0.2)  # Adjust the vertical position of the labels closer to the colorbar

    # plot horizantal bars
    if "Z_" == score[:2]:
        df_copy = df.copy()
        for targ_list in targets_to_choose_best:
            targ_index = [targ for targ in targ_list if targ in pivot_all.columns]
            best_scores = pivot_all[targ_index].max(1)
            pivot_all[targ_list[0]] = best_scores
            pivot_all[targ_list[1:]] = 0 
        pivot_all["sum z>0"] = pivot_all.where(pivot_all > 0).sum(1)
        plt.barh(g.get_yticks(),width=pivot_all["sum z>0"],color=h_bar_color)
        ax_yDist.set_xlabel("Sum Z>0")
    else:
        cbar.set_label("Sum",size=ticklabel_size*1.2,
            ha='right',y=-0.1, rotation=0, va='bottom')
        for targ_list in targets_to_choose_best:
            targ_index = [targ for targ in targ_list if targ in pivot_all.columns]
            best_scores = pivot_all[targ_index].max(1)
            pivot_all[targ_list[0]] = best_scores
            pivot_all[targ_list[1:]] = 0 
        pivot_all["Sum"] = pivot_all.sum(1)
        plt.barh(g.get_yticks(),width=pivot_all["Sum"],color=h_bar_color)
        ax_yDist.set_xlabel("Sum")
    # TODO if want to plot another
    #zrna_scores = get_group_score(df_copy)
    #plt.barh(g.get_yticks(),width=zrna_scores.Z_rna,edgecolor="orange", fill=False)
    t=ax_yDist.set_yticklabels(labels=g.get_yticklabels(),rotation = 0,size=y_text_size)
    ax_yDist.spines['right'].set_visible(False)
    ax_yDist.spines['top'].set_visible(False)
    
    #print(pivot_all["sum z>0"])
    #from matplotlib.lines import Line2D
    #custom_lines = [Line2D([0], [0], color="green", lw=4),
                    #Line2D([0], [0], color="orange", lw=4)]

    #ax_yDist.legend(custom_lines, ['Z EM', 'Z RNA'],loc='lower right')

    if 'AF3-server - 304' in pivot_all.index:
        # Plot AF3 baseline
        x_coord = pivot_all["sum z>0" if "Z_" in score else "Sum"]['AF3-server - 304'] 
        ax_yDist.axvline(x=x_coord, color='gold', linestyle='--', linewidth=1, label='AF3-server - 304')

    ax_main_pos = ax_main.get_position()  # Get position of the main axis
    ax_yDist.set_position([ax_yDist.get_position().x0, ax_main_pos.y0, ax_yDist.get_position().width, ax_main_pos.height])

    if savefig is not None:
        plt.savefig(f'{savefig}.png',dpi=400, bbox_inches='tight')
        plt.savefig(f'{savefig}.svg',dpi=400, bbox_inches='tight')

    return rank_order



def plot_barplot_SS(preds_top1, scores, palette,figsize=(10, 6),top_N=None,y_highlight=[],save=None,
                ticklabel_size=5,title='',yrot=0,rank_type='sum',targets_participating=None,
                group_col='group',targets_to_choose_best=TARGETS_CHOOSE_BEST):
    df_copy = preds_top1.copy()

    if rank_type == 'sum':
        group_score = []
        for score in scores:
            group_score.append(get_group_score(df_copy, agg="sum", score=score,
                    targets_to_choose_best=targets_to_choose_best,group_col=group_col,
                   fill_nan=1,fill_not_present=0).sort_values(score,ascending=False).set_index(group_col))
        group_score = pd.concat(group_score)
        print(group_score)
    elif rank_type == 'mean':
        group_score = get_group_score(df_copy, agg="mean", score=score, # mean ignores nans so if 6 values, 2 nan, will noramlize sum by 4
                    targets_to_choose_best=targets_to_choose_best,group_col=group_col,
                   fill_nan=np.nan,fill_not_present=0).sort_values(score,ascending=False)
    rank_order = group_score.index.values

    if rank_type == 'sum':
        group_sums = df_copy.groupby('group')[scores].sum().sum(axis=1).sort_values(ascending=False)
        rank_order = group_sums.index.values
        grouped = df_copy.groupby('group')[scores].sum().reset_index().copy()
        
        # for each score need to normalize by number of targets in category
        for score in scores:
            grouped[score] = grouped[score]/len(targets_participating[score])
    elif rank_type=='sum_norm':
        # rank by sum, with nan being a score of 1
        rank_order_sum_nan_one = (df_copy.groupby('group')
            .apply(lambda x: x[scores].fillna(1).sum()).sum(axis=1)  # Replace NaN with 1 and then sum
            .reset_index(name='total_score').sort_values('total_score', ascending=False))
        rank_order = rank_order_sum_nan_one['group'].values
        grouped = (df_copy.groupby('group')
            .apply(lambda x: x[scores].fillna(1).sum())).reset_index().copy()
        # for each score need to min max normalize. Max is total targets, min is total targets - targets part.
        for score in scores:
            max_score = len(df_copy.target.unique())
            min_score = max_score - len(targets_participating[score])
            grouped[score] = ( grouped[score]-min_score) / (max_score-min_score)
    elif rank_type=='mean':
        # rank by nan mean (aka sum and then normalie by number of non nan values)
        rank_order_nan_mean = (df_copy.groupby('group')
            .apply(lambda x: x[scores].mean(skipna=True)).sum(axis=1)
            .reset_index(name='mean_score').sort_values('mean_score', ascending=False)) 
        rank_order = rank_order_nan_mean['group'].values
        grouped = (df_copy.groupby('group')
            .apply(lambda x: x[scores].mean(skipna=True))).reset_index().copy()

    if not top_N is None:
        rank_order = rank_order[:top_N]
    preds_long = grouped.melt(id_vars=['group'], value_vars=scores, 
        var_name='score_type', value_name='score')
    

    fig, ax = plt.subplots(figsize=figsize)  
    sns.barplot(preds_long,x='score',y='group',order=rank_order,
                hue='score_type',hue_order=scores,palette=palette,legend=False,
                width=0.8,linewidth=0.5, edgecolor="black")
    ax.spines[['top','right']].set_visible(False)
    #ax.set_xlim(0,1)
    ax.set_title(title)
    t=ax.set_yticklabels(labels=ax.get_yticklabels(),rotation = yrot,size=ticklabel_size)
    t=ax.set_xticklabels(labels=ax.get_xticklabels(),rotation = 0,size=ticklabel_size+1)
    ax.set_ylabel("Group")
    ax.set_xlabel("Sum F1-score")
    for label in ax.get_yticklabels():
        if label.get_text() in y_highlight:
            label.set_color('grey')
    if save:
        plt.savefig(f"{save}.svg",dpi=400, bbox_inches='tight', transparent=True)
        plt.savefig(f"{save}.png",dpi=400, bbox_inches='tight', transparent=True)

'''