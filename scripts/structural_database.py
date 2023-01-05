# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 09:39:14 2022

@author: DELL
"""


import os
import pickle
import numpy as np
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem

source_path = "E:/Software/MSFINDER ver 3.52/MS-FINDER ver. 3.52 Windows/Resources/MsfinderStructureDB-VS15.esd"
save_path = "D:/DeepMASS2_GUI/data/DeepMassStructureDB-v1.0.csv"
mf_database = pd.read_csv(source_path, sep='\t')
path_data = 'D:/All_MSDatabase'

# Add GNPS
outfile = os.path.join(path_data, 'GNPS_all/ALL_GNPS_220601_positive_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums = pickle.load(file)

outfile = os.path.join(path_data, 'GNPS_all/ALL_GNPS_220601_negative_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums += pickle.load(file)

database, short_keys = [], []
for i, s in enumerate(tqdm(spectrums)):
    short_key = s.metadata['inchikey'][:14]
    if short_key == '':
        continue
    if (short_key not in mf_database['Short InChIKey'].values) and (short_key not in short_keys):
        smi = s.metadata['smiles']
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        formula = AllChem.CalcMolFormula(mol)
        exactmass = AllChem.CalcExactMolWt(mol)
        try:
            database.append([s.metadata['compound_name'],
                             s.metadata['smiles'],
                             formula,
                             s.metadata['inchikey'],
                             short_key,
                             exactmass,
                             'GNPS{}'.format(i)])
            short_keys.append(short_key)
        except:
            pass

database = pd.DataFrame(database, columns=['Title', 'SMILES', 'Formula', 'InChIkey', 'Short InChIKey', 'Exact mass', 'GNPS'])
mf_database = mf_database.append(database)


# Add NIST
outfile = os.path.join(path_data, 'NIST2020/ALL_NIST20_positive_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums = pickle.load(file)

outfile = os.path.join(path_data, 'NIST2020/ALL_NIST20_negative_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums += pickle.load(file)

database, short_keys = [], []
for i, s in enumerate(tqdm(spectrums)):
    short_key = s.metadata['inchikey'][:14]
    if short_key == '':
        continue
    if (short_key not in mf_database['Short InChIKey'].values) and (short_key not in short_keys):
        smi = s.metadata['smiles']
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        formula = AllChem.CalcMolFormula(mol)
        exactmass = AllChem.CalcExactMolWt(mol)
        try:
            database.append([s.get('synon', None),
                             s.get('smiles', None),
                             formula,
                             s.get('inchikey', None),
                             short_key,
                             exactmass,
                             s.get('nistno', None)])
            short_keys.append(short_key)
        except:
            pass
database = pd.DataFrame(database, columns=['Title', 'SMILES', 'Formula', 'InChIkey', 'Short InChIKey', 'Exact mass', 'NIST'])
mf_database = mf_database.append(database)

mf_database = mf_database.sort_values(by = 'Exact mass')
mf_database.to_csv(save_path, index = False)
