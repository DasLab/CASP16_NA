import json
from glob import glob
import numpy as np
import pandas as pd


def interface_list_to_identity(interfaces,name):
    proteins = [chr(x) for x in range(65,90)]
    RNA = [str(x) for x in range(15)] # and DNA....
    # O to Q DNA? M1268
    # M1297_v1 b in RNA
    if 'M1268' in name or 'M0268' in name:
        proteins = [chr(x) for x in range(65,79)]+[chr(x) for x in range(82,90)]
        RNA = [chr(x) for x in range(79,82)]
    if 'M1297' in name:
        RNA = ['b']
        proteins = ['c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k']

    interface_types = []
    for interface in interfaces:
        if interface[0] not in proteins+RNA or interface[1] not in proteins+RNA:
            print("ERROR",interfaces)
        if interface[0] in proteins and interface[1] in proteins:
            interface_types.append('protein-protein')
        elif interface[0] in RNA and interface[1] in RNA:
            interface_types.append('NA-NA')
        else:
            interface_types.append('NA-protein')
    return interface_types

def get_interface_residues(interfaces,residues):
    int_dict = {tuple(x):[] for x in interfaces}
    for res in residues:
        int_dict[(res[0].split('.')[0],res[1].split('.')[0])].extend(res)
    res_list = []
    for interface in interfaces:
        res_list.append(list(set(int_dict[tuple(interface)])))
    return res_list

def get_interface_lddt_size(interfaces,residues,lddt_scores,chain_mapping):
    res_list = get_interface_residues(interfaces,residues)
    lddt = []
    #print(lddt_scores)
    for ress in res_list:
        lddts = []
        for res in ress:
            if res[0] not in chain_mapping:
                lddts.append(0)
            elif chain_mapping[res[0]]+res[1:] not in lddt_scores:
                lddts.append(0)
            elif lddt_scores[chain_mapping[res[0]]+res[1:]] is None:
                lddts.append(0)
            else:
                lddts.append(lddt_scores[chain_mapping[res[0]]+res[1:]])
        lddts[lddts == None] = np.nan
        lddts = np.nan_to_num(lddts)
        lddt.append(sum(lddts)/len(lddts))
    return np.array(lddt), np.array([len(x) for x in res_list])

def normalize_score(score,size,method='uniform'):
    if method == 'uniform':
        return (score*size).sum() / size.sum()
    elif method == 'log10':
        return (score*np.log10(size/2)).sum() / np.log10(size/2).sum()
    else:
        print(f'method {method} not implemented')
        return None
    

results = []
for file in sorted(glob("ost/M*/*.json")):
    with open(file, "r") as f:
        # print(file)
        data = json.load(f)
        
        # ignore anyone without all chains
        if data['status']=='FAILURE':
            continue
        if len(data['reference_chains']) > len(data['model_chains']):
            continue

        model = file.split('/')[-1][:-5]
        target = model.split('TS')[0]
        iterface_identity = interface_list_to_identity(data['contact_reference_interfaces'],file)
        interface_ics = np.array(data['per_interface_ics'])
        interface_ips = np.array(data['per_interface_ips'])
        lddt,size = get_interface_lddt_size(data['contact_reference_interfaces'],data['reference_contacts'],data['local_lddt'],data['chain_mapping'])

        interface_ics[interface_ics == None] = np.nan
        interface_ips[interface_ips == None] = np.nan
        lddt[lddt == None] = np.nan
        
        interface_ics[interface_ics!=interface_ics] = 0
        interface_ips[interface_ips!=interface_ips] = 0
        lddt[lddt!=lddt] = 0
        
        na_na = np.array(iterface_identity) == 'NA-NA'
        pro_pro = np.array(iterface_identity) == 'protein-protein'
        na_pro = np.array(iterface_identity) == 'NA-protein'
        na_containing = na_na | na_pro

        new_data = ['all',model,target]
        for score in [interface_ics,interface_ips,lddt]:
            new_data.append (score.mean())
            for method in ['uniform','log10']:
                new_data.append(normalize_score(score,size,method=method))
            
        results.append(new_data)
        for name,select in zip(['NA-NA','NA-protein','protein-protein','NA-containing'],[na_na,na_pro,pro_pro,na_containing]):
            new_data = [name,model,target]
            for score in [interface_ics,interface_ips,lddt]:
                new_data.append (score.mean())
                for method in ['uniform','log10']:
                    new_data.append(normalize_score(score,size,method=method))
            results.append(new_data)

df = pd.DataFrame(results,columns=['interface_type','model','target','ics_naive_mean','ics_mean','ics_log10_mean','ips_naive_mean','ips_mean','ips_log10_mean','ilddt_naive_mean','ilddt_mean','ilddt_log10_mean',])

df.to_csv('casp16_domain_score_rck.csv',index=False)