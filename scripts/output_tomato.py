# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 08:25:08 2023

@author: DELL
"""

import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from matchms.importing import load_from_mgf

dir_pos = 'example/Tomato/identify/positive'
msdial_pos = pd.read_csv('example/Tomato/tomato_positive_msdial.csv')

spectrums_pos = [s for s in load_from_mgf('example/Tomato/ms_ms_tomato_all_positive.mgf')]
index_pos = [[s.get('compound_name'), s.get('peak_index')] for s in spectrums_pos]
index_pos = pd.DataFrame(index_pos, columns =['index', 'peak'])

res_pos = os.listdir(dir_pos)
keep_pos, identify_pos, inchikey_pos, smiles_pos = [], [], [], []
for f in tqdm(res_pos):
    peak = f.split('.')[0]
    w = np.where(index_pos['index'] == peak)[0]
    if len(w) == 0:
        continue
    else:
        w = w[0]
    ind = index_pos.loc[w, 'peak']
    f1 = dir_pos + '/{}'.format(f)
    ann = pd.read_csv(f1)
    if len(ann) == 0:
        continue

    k = int(ind.split('_')[1])
    keep_pos.append(k)
    
    if 'Title' in ann.columns:
        name = ann.loc[0, 'Title']
        inchikey_pos.append(ann.loc[0, 'InChIKey'])
        smiles_pos.append(ann.loc[0, 'CanonicalSMILES'])
        identify_pos.append(name)
    else:
        name = 'PubChem_{}'.format(ann.loc[0, 'CID'])
        inchikey_pos.append(ann.loc[0, 'InChIKey'])
        smiles_pos.append(ann.loc[0, 'CanonicalSMILES'])
        identify_pos.append(name)

keep_pos = np.array(keep_pos)
output_pos = msdial_pos.loc[keep_pos, :]
output_pos['Metabolite name'] = identify_pos
output_pos['INCHIKEY'] = inchikey_pos
output_pos['SMILES'] = smiles_pos
output_pos.to_csv('example/Tomato/tomato_positive_output.csv')
    

dir_neg = 'example/Tomato/identify/negative'
msdial_neg = pd.read_csv('example/Tomato/tomato_negative_msdial.csv')

spectrums_neg = [s for s in load_from_mgf('example/Tomato/ms_ms_tomato_all_negative.mgf')]
index_neg = [[s.get('compound_name'), s.get('peak_index')] for s in spectrums_neg]
index_neg = pd.DataFrame(index_neg, columns =['index', 'peak'])

res_neg = os.listdir(dir_neg)
keep_neg, identify_neg, inchikey_neg, smiles_neg = [], [], [], []
for f in tqdm(res_neg):
    peak = f.split('.')[0]
    w = np.where(index_neg['index'] == peak)[0]
    if len(w) == 0:
        continue
    else:
        w = w[0]
    ind = index_neg.loc[w, 'peak']
    f1 = dir_neg + '/{}'.format(f)
    ann = pd.read_csv(f1)
    if len(ann) == 0:
        continue

    k = int(ind.split('_')[1])
    keep_neg.append(k)
    
    if 'Title' in ann.columns:
        name = ann.loc[0, 'Title']
        inchikey_neg.append(ann.loc[0, 'InChIKey'])
        smiles_neg.append(ann.loc[0, 'CanonicalSMILES'])
        identify_neg.append(name)
    else:
        name = 'PubChem_{}'.format(ann.loc[0, 'CID'])
        inchikey_neg.append(ann.loc[0, 'InChIKey'])
        smiles_neg.append(ann.loc[0, 'CanonicalSMILES'])
        identify_neg.append(name)

keep_neg = np.array(keep_neg)
output_neg = msdial_neg.loc[keep_neg, :]
output_neg['Metabolite name'] = identify_neg
output_neg['INCHIKEY'] = inchikey_neg
output_neg['SMILES'] = smiles_neg
output_neg.to_csv('example/Tomato/tomato_negative_output.csv')

