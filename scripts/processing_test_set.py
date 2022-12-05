# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 09:59:18 2022

@author: DELL
"""

import os
import numpy as np
import pickle

from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem

from matchms.exporting import save_as_mgf, save_as_msp

inchikey_test = np.load('data/inchikey_test.npy', allow_pickle=True)

path_data = 'D:/All_MSDatabase'
outfile = os.path.join(path_data, 'GNPS_all/ALL_GNPS_220601_positive_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums = pickle.load(file)

outfile = os.path.join(path_data, 'In_House/ALL_Inhouse_positive_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums += pickle.load(file)

outfile = os.path.join(path_data, 'NIST2020/ALL_NIST20_positive_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums += pickle.load(file)

outfile = os.path.join(path_data, 'GNPS_all/ALL_GNPS_220601_negative_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums += pickle.load(file)

outfile = os.path.join(path_data, 'In_House/ALL_Inhouse_negative_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums += pickle.load(file)

outfile = os.path.join(path_data, 'NIST2020/ALL_NIST20_negative_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums += pickle.load(file)

spectrums = [s for s in tqdm(spectrums) if s.get('inchikey')[:14] in inchikey_test]


testing_spectrums = []
for s in tqdm(spectrums):
    if len(s.mz[s.mz <= 1000]) < 5:
        continue
    if 'smiles' not in list(s.metadata.keys()):
        continue
    if s.metadata['smiles'] == '':
        continue
    try:
        mol = Chem.MolFromSmiles(s.metadata['smiles'])
        wt = AllChem.CalcExactMolWt(mol)
        smi = Chem.MolToSmiles(mol, isomericSmiles=False)
        s = s.set('compound_name', 'challenge_{}'.format(len(testing_spectrums)))
        s = s.set('smiles', smi)
        s = s.set('parent_mass', wt)
    except:
        continue
    if 'ionmode' not in list(s.metadata.keys()):
        continue  

    keys = [s for s in s.metadata.keys() if s in ['compound_name', 'smiles', 'inchikey', 'precursor_mz', 'adduct', 'parent_mass', 'ionmode', 'charge']]
    s.metadata = {k: s.metadata[k] for k in keys}
    testing_spectrums.append(s)


# save as mgf
save_as_mgf(list(testing_spectrums), 'example/Test/all_test_set.mgf')

# save as msp individually
for i, s in enumerate(testing_spectrums):
    path = 'example/Test/msp/challenge_{}.msp'.format(i)
    save_as_msp([testing_spectrums[i]], path)
    
    with open(path, encoding = 'utf-8') as msp:
        lines = msp.readlines()
        lines = [l.replace('_', '') for l in lines]
        lines = [l.replace('ADDUCT', 'PRECURSORTYPE') for l in lines]
    with open(path, 'w') as msp:
        msp.writelines(lines)

    # only for ms-finder
    path_msfinder = path.replace('/msp/', '/msfinder/')
        
    # exclude large compound, as ms-finder is very slow for them
    if s.metadata['parent_mass'] >= 850:
        continue
    with open(path_msfinder, 'w') as msp:
        msp.writelines(lines)