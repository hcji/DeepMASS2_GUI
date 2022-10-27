# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 08:46:57 2022

@author: DELL
"""

import os
import numpy as np
import pickle
import hnswlib
from tqdm import tqdm
from rdkit import Chem

from ms2deepscore import MS2DeepScore
from ms2deepscore.models import load_model

path_data = 'D:/All_MSDatabase'

outfile = os.path.join(path_data, 'GNPS_all/ALL_GNPS_220601_positive_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums = pickle.load(file)

outfile = os.path.join(path_data, 'In_House/ALL_In_House_positive_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums += pickle.load(file)

outfile = os.path.join(path_data, 'NIST2017/ALL_NIST17_positive_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums += pickle.load(file)


reference = []
for s in tqdm(spectrums):
    if len(s.mz[s.mz <= 1000]) < 5:
        continue
    if 'smiles' not in list(s.metadata.keys()):
        continue
    try:
        smi = Chem.MolToSmiles(Chem.MolFromSmiles(s.metadata['smiles']), isomericSmiles=False)
        s.metadata['smiles'] = smi
    except:
        continue
    if 'ionmode' not in list(s.metadata.keys()):
        continue  
    if s.metadata['ionmode'] == 'positive':
        reference.append(s)


model = load_model(os.path.join(path_data, 'MS2DeepScore_allGNPSpositive.hdf5'))
ms2deepscore = MS2DeepScore(model)
reference_vector = ms2deepscore.calculate_vectors(reference)

xb = np.array(reference_vector).astype('float32')
xb_len =  np.linalg.norm(xb, axis=1, keepdims=True)
xb = xb/xb_len
dim = 200
num_elements = len(xb)
ids = np.arange(num_elements)

p = hnswlib.Index(space = 'l2', dim = dim)
p.init_index(max_elements = num_elements, ef_construction = 800, M = 64)
p.add_items(xb, ids)
p.set_ef(300)
p.save_index('data/references_index_positive.bin')
np.save('data/references_spectrums_positive.npy', reference)
