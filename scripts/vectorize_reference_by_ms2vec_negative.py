# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 08:46:57 2022

@author: DELL
"""

import os
import random
import numpy as np
import pickle
import hnswlib
import gensim
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem

from matchms.importing import load_from_mgf
from spec2vec import SpectrumDocument
from spec2vec.vector_operations import calc_vector

inchikey_test = np.load('data/inchikey_test.npy', allow_pickle=True)
spectrums_casmi = [s for s in load_from_mgf('example/CASMI/all_casmi.mgf')]
inchikey_casmi = list(set([s.get('inchikey')[:14] for s in spectrums_casmi]))
inchikey_test = [s for s in inchikey_test if s not in inchikey_casmi]

# negative
path_data = 'D:/All_MSDatabase'

outfile = os.path.join(path_data, 'GNPS_all/ALL_GNPS_220601_negative_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums = pickle.load(file)

outfile = os.path.join(path_data, 'In_House/ALL_Inhouse_negative_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums += pickle.load(file)

outfile = os.path.join(path_data, 'NIST2020/ALL_NIST20_negative_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums += pickle.load(file)

spectrums = [s for s in tqdm(spectrums)]
random.shuffle(spectrums)


reference, test = [], []
count_test = np.zeros(len(inchikey_test))
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
        formula = AllChem.CalcMolFormula(mol)
        s = s.set('formula', formula)
        s = s.set('smiles', smi)
        s = s.set('parent_mass', wt)
    except:
        continue
    if 'ionmode' not in list(s.metadata.keys()):
        continue
    if s.metadata['ionmode'] == 'negative':
        keys = [s for s in s.metadata.keys() if s in ['compound_name', 'formula', 'smiles', 'inchikey', 'precursor_mz', 'adduct', 'parent_mass', 'ionmode', 'charge']]
        s.metadata = {k: s.metadata[k] for k in keys}
        inchikey = s.get('inchikey')[:14]
        if inchikey in inchikey_casmi:
            continue
        if inchikey in inchikey_test:
            i = inchikey_test.index(inchikey)
            if count_test[i] == 0:
                count_test[i] += 1
                test.append(s)
                continue
        reference.append(s)


file = 'model/Ms2Vec_allGNPSnegative.hdf5_iter_30.model'
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
p.save_index('data/references_index_negative_spec2vec.bin')
pickle.dump(reference, open('data/references_spectrums_negative.pickle', "wb"))
pickle.dump(test, open('example/Test/test_spectrums_negative.pickle', "wb"))