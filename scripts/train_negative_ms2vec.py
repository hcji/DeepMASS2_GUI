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

path_data = 'D:/All_MSDatabase'

outfile = os.path.join(path_data, 'GNPS_all/ALL_GNPS_220601_negative_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums = pickle.load(file)

outfile = os.path.join(path_data, 'In_House/ALL_In_House_negative_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums += pickle.load(file)

outfile = os.path.join(path_data, 'NIST2017/ALL_NIST17_negative_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums += pickle.load(file)

print("number of spectra:", len(spectrums))


documents = [SpectrumDocument(s, n_decimals=2) for s in spectrums]

model_file = "model/Ms2Vec_allGNPSnegative.hdf5"
iterations = [30, 40, 50]
model = train_new_word2vec_model(documents, iterations=iterations, filename=model_file,
                                 workers=8, progress_logger=True)