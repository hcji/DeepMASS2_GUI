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
save_path = "D:/DeepMASS2_GUI/data/MsfinderStructureDB-VS15-plus-GNPS.csv"

mf_database = pd.read_csv(source_path, sep='\t')

spectrums = list(np.load('data/references_spectrums_positive.npy', allow_pickle=True))
spectrums += list(np.load('data/references_spectrums_negative.npy', allow_pickle=True))

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
mf_database = mf_database.sort_values(by = 'Exact mass')
mf_database.to_csv(save_path, index = False)
