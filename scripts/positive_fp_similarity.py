# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 14:11:47 2022

@author: DELL
"""

import os
import numpy as np
import pickle
from tqdm import tqdm

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

print("number of spectra:", len(spectrums))



inchikeys_list = []
for s in spectrums:
    inchikeys_list.append(s.get("inchikey"))

inchikeys14_array = np.array([x[:14] for x in inchikeys_list])
inchikeys14_unique = list({x[:14] for x in inchikeys_list})

from collections import Counter 
  
def most_frequent(List): 
    occurence_count = Counter(List) 
    return occurence_count.most_common(1)[0][0] 

inchi_list = []
for s in spectrums:
    inchi_list.append(s.get("inchi"))

inchi_array = np.array(inchi_list)


inchi_mapping = []
ID_mapping = []

for inchikey14 in inchikeys14_unique:
    idx = np.where(inchikeys14_array == inchikey14)[0]
    
    inchi = most_frequent([spectrums[i].get("inchi") for i in idx])
    inchi_mapping.append(inchi)
    ID = idx[np.where(inchi_array[idx] == inchi)[0][0]]
    ID_mapping.append(ID)


import pandas as pd
metadata = pd.DataFrame(list(zip(inchikeys14_unique, inchi_mapping, ID_mapping)), columns=["inchikey", "inchi", "ID"])

from matchms.filtering.add_fingerprint import add_fingerprint

for i in tqdm(metadata.ID.values):
    spectrums[i] = add_fingerprint(spectrums[i], fingerprint_type="daylight", nbits=2048)


from matchms.similarity import FingerprintSimilarity

spectrums_represent = [spectrums[i] for i in metadata.ID.values]
similarity_measure = FingerprintSimilarity(similarity_measure="jaccard")
scores_mol_similarity = similarity_measure.matrix(spectrums_represent, spectrums_represent)

tanimoto_df = pd.DataFrame(scores_mol_similarity, columns=metadata.inchikey.values, index=metadata.inchikey.values)


filename = os.path.join(path_data, "ALL_positive_tanimoto_scores.pickle")
tanimoto_df.to_pickle(filename)
