# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 14:28:08 2022

@author: DELL
"""

import os
import numpy as np
import pickle
from tqdm import tqdm
from spec2vec import SpectrumDocument
from spec2vec.model_building import train_new_word2vec_model

inchikey_test = np.load('data/inchikey_test.npy', allow_pickle=True)

# positive
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

spectrums = [s for s in tqdm(spectrums) if s.get('inchikey')[:14] not in inchikey_test]

documents = [SpectrumDocument(s, n_decimals=2) for s in spectrums]
model_file = "model/Ms2Vec_allGNPSpositive.hdf5"
iterations = [30, 40]
model = train_new_word2vec_model(documents, iterations=iterations, filename=model_file,
                                 workers=8, progress_logger=True)



# negative
outfile = os.path.join(path_data, 'GNPS_all/ALL_GNPS_220601_negative_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums = pickle.load(file)

outfile = os.path.join(path_data, 'In_House/ALL_Inhouse_negative_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums += pickle.load(file)

outfile = os.path.join(path_data, 'NIST2020/ALL_NIST20_negative_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums += pickle.load(file)

spectrums = [s for s in tqdm(spectrums) if s.get('inchikey')[:14] not in inchikey_test]

documents = [SpectrumDocument(s, n_decimals=2) for s in spectrums]
model_file = "model/Ms2Vec_allGNPSnegative.hdf5"
iterations = [30, 40]
model = train_new_word2vec_model(documents, iterations=iterations, filename=model_file,
                                 workers=8, progress_logger=True)


