# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 08:46:57 2022

@author: DELL
"""

import os
import numpy as np
import pickle
import hnswlib
import gensim
from tqdm import tqdm
from rdkit import Chem

from spec2vec import SpectrumDocument
from spec2vec.vector_operations import calc_vector

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
    if s.metadata['smiles'] == '':
        continue
    try:
        smi = Chem.MolToSmiles(Chem.MolFromSmiles(s.metadata['smiles']), isomericSmiles=False)
        s.metadata['smiles'] = smi
    except:
        continue
    if 'ionmode' not in list(s.metadata.keys()):
        continue  
    if s.metadata['ionmode'] == 'positive':
        keys = [s for s in s.metadata.keys() if s in ['compound_name', 'smiles', 'inchikey', 'precursor_mz', 'adduct', 'parent_mass', 'ionmode', 'charge']]
        s.metadata = {k: s.metadata[k] for k in keys}
        reference.append(s)


file = 'model/Ms2Vec_allGNPSpositive.hdf5_iter_30.model'
model = gensim.models.Word2Vec.load(file)
calc_ms2vec_vector = lambda x: calc_vector(model, SpectrumDocument(x, n_decimals=2))

reference_vector = []
for s in tqdm(reference):
    reference_vector.append(calc_ms2vec_vector(s))

xb = np.array(reference_vector).astype('float32')
xb_len =  np.linalg.norm(xb, axis=1, keepdims=True)
xb = xb/xb_len
dim = 300
num_elements = len(xb)
ids = np.arange(num_elements)

p = hnswlib.Index(space = 'l2', dim = dim)
p.init_index(max_elements = num_elements, ef_construction = 800, M = 64)
p.add_items(xb, ids)
p.set_ef(300)
p.save_index('data/references_index_positive_spec2vec.bin')