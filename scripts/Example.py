# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 08:49:29 2023

@author: DELL
"""


import umap
import hnswlib
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import scale
from sklearn.decomposition import PCA
from matchms.importing import load_from_mgf
from gensim.models import Word2Vec
from rdkit import Chem
from rdkit.Chem import DataStructs, AllChem, Draw
from spec2vec import SpectrumDocument
from spec2vec.vector_operations import calc_vector

from core.identification import identify_unknown, match_spectrum

database = pd.read_csv('data/DeepMassStructureDB-v1.0.csv')
spectrums = [s for s in load_from_mgf("D:/DeepMASS2_Data_Processing/Example/CASMI/all_casmi.mgf")]

# Example 1
s = spectrums[255]
print(s.metadata)

model = Word2Vec.load("model/Ms2Vec_allGNPSnegative.hdf5")
p = hnswlib.Index(space='l2', dim=300) 
p.load_index('data/references_index_negative_spec2vec.bin')
with open('data/references_spectrums_negative.pickle', 'rb') as file:
    references = pickle.load(file)
references = np.array(references)
precursors = [s.get('precursor_mz') for s in references]
precursors = np.array(precursors)

s_metadata = s.metadata
s_matchms = match_spectrum(s, precursors, references)
s_matchms_metadata = s_matchms.metadata
s_deepmass = identify_unknown(s, p, model, references, database)
s_deepmass_metadata = s_deepmass.metadata


get_mol_fingerprint = lambda x: AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(x), radius=2)
get_mol_similarity = lambda x, y: DataStructs.FingerprintSimilarity(x, y)
calc_ms2vec_vector = lambda x: calc_vector(model, SpectrumDocument(x, n_decimals=2))
deepmass_candidate = s_deepmass_metadata['annotation']['CanonicalSMILES']
deepmass_reference = s_deepmass_metadata['reference']

deepmass_candidate_vector = []
for r in deepmass_candidate:
    # a = calc_ms2vec_vector(s)
    b = get_mol_fingerprint(r)
    deepmass_candidate_vector.append(list(b))
deepmass_candidate_vector = np.array(deepmass_candidate_vector)

deepmass_reference_vector = []
for r in deepmass_reference:
    # a = calc_ms2vec_vector(r)
    b = get_mol_fingerprint(r.get('smiles'))
    deepmass_reference_vector.append(list(b))
deepmass_reference_vector = np.array(deepmass_reference_vector)

X = np.vstack((deepmass_candidate_vector, deepmass_reference_vector))
pca = umap.UMAP(n_components = 2)
X_r = pca.fit_transform(scale(X))

Draw.MolToFile(Chem.MolFromSmiles(s.get('smiles')), 'temp/temp1.png')

a, b = len(deepmass_candidate_vector), len(deepmass_reference_vector)
plt.figure(dpi = 300, figsize = (3.5, 3.5))
plt.scatter(X_r[1:a,0], X_r[1:a,1], color = 'green', alpha = 0.5, label = 'Candidates')
plt.scatter(X_r[a:a+10,0], X_r[a:a+10,1], color = 'blue', alpha = 0.5, label = 'Top 10 Reference')
plt.scatter(X_r[0,0], X_r[0,1], color = 'red', alpha = 0.8, label = 'True Annotation')
plt.scatter(X_r[a+9:a+10,0], X_r[a+9:a+10,1], color = 'orange', alpha = 0.5, label = 'Top 10 Reference')
plt.xlabel('Dim 1')
plt.ylabel('Dim 2')
plt.xlim(40, 50)
plt.ylim(20, 25)
plt.legend()



# Example 2
s = spectrums[368]
print(s.metadata)

model = Word2Vec.load("model/Ms2Vec_allGNPSpositive.hdf5")
p = hnswlib.Index(space='l2', dim=300) 
p.load_index('data/references_index_positive_spec2vec.bin')
with open('data/references_spectrums_positive.pickle', 'rb') as file:
    references = pickle.load(file)
references = np.array(references)
precursors = [s.get('precursor_mz') for s in references]
precursors = np.array(precursors)

model = Word2Vec.load("model/Ms2Vec_allGNPSpositive.hdf5")
p = hnswlib.Index(space='l2', dim=300) 
p.load_index('data/references_index_positive_spec2vec.bin')
with open('data/references_spectrums_positive.pickle', 'rb') as file:
    references = pickle.load(file)
references = np.array(references)
precursors = [s.get('precursor_mz') for s in references]
precursors = np.array(precursors)

s_metadata = s.metadata
s_matchms = match_spectrum(s, precursors, references)
s_matchms_metadata = s_matchms.metadata
s_deepmass = identify_unknown(s, p, model, references, database)
s_deepmass_metadata = s_deepmass.metadata


get_mol_fingerprint = lambda x: AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(x), radius=2)
get_mol_similarity = lambda x, y: DataStructs.FingerprintSimilarity(x, y)
calc_ms2vec_vector = lambda x: calc_vector(model, SpectrumDocument(x, n_decimals=2))
deepmass_candidate = s_deepmass_metadata['annotation']['CanonicalSMILES']
deepmass_reference = s_deepmass_metadata['reference']

deepmass_candidate_vector = []
for r in deepmass_candidate:
    b = get_mol_fingerprint(r)
    deepmass_candidate_vector.append(list(b))
deepmass_candidate_vector = np.array(deepmass_candidate_vector)

deepmass_reference_vector = []
for r in deepmass_reference:
    b = get_mol_fingerprint(r.get('smiles'))
    deepmass_reference_vector.append(list(b))
deepmass_reference_vector = np.array(deepmass_reference_vector)

X = np.vstack((deepmass_candidate_vector, deepmass_reference_vector))
pca = umap.UMAP(n_components = 2)
X_r = pca.fit_transform(scale(X))

Draw.MolToFile(Chem.MolFromSmiles(s.get('smiles')), 'temp/temp.png')

a, b = len(deepmass_candidate_vector), len(deepmass_reference_vector)
plt.figure(dpi = 300, figsize = (3.5, 3.5))
plt.scatter(X_r[1:a,0], X_r[1:a,1], color = 'green', alpha = 0.5, label = 'Candidates')
plt.scatter(X_r[0,0], X_r[0,1], color = 'red', alpha = 0.8, label = 'True Annotation')
plt.scatter(X_r[a:a+10,0], X_r[a:a+10,1], color = 'blue', alpha = 0.5, label = 'Top 10 Reference')
plt.scatter(X_r[a+9:a+10,0], X_r[a+9:a+10,1], color = 'orange', alpha = 0.5, label = 'Top 10 Reference')
plt.xlabel('Dim 1')
plt.ylabel('Dim 2')
plt.legend()