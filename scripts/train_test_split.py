# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 08:10:18 2022

@author: DELL
"""

import os
import numpy as np
import pickle
import random
from tqdm import tqdm
from matchms.importing import load_from_mgf

# casmi
spectrums = [s for s in load_from_mgf('example/CASMI/all_casmi.mgf')]
inchikey_casmi = list(set([s.get('inchikey')[:14] for s in spectrums]))


# positive
path_data = 'D:/All_MSDatabase'

outfile = os.path.join(path_data, 'GNPS_all/ALL_GNPS_220601_positive_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums = pickle.load(file)

outfile = os.path.join(path_data, 'In_House/ALL_In_House_positive_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums += pickle.load(file)

outfile = os.path.join(path_data, 'NIST2020/ALL_NIST20_positive_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums += pickle.load(file)

print("number of spectra:", len(spectrums))
inchikey_positive = set([s.get('inchikey')[:14] for s in spectrums])
inchikey_positive = random.sample(inchikey_positive, k = 1500)


# negative
outfile = os.path.join(path_data, 'GNPS_all/ALL_GNPS_220601_negative_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums = pickle.load(file)

outfile = os.path.join(path_data, 'In_House/ALL_In_House_negative_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums += pickle.load(file)

outfile = os.path.join(path_data, 'NIST2020/ALL_NIST20_negative_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums += pickle.load(file)

print("number of spectra:", len(spectrums))
inchikey_negative = set([s.get('inchikey')[:14] for s in spectrums])
inchikey_negative = random.sample(inchikey_negative, k = 1500)

inchikey_test = list(set(inchikey_casmi + inchikey_positive + inchikey_negative))
np.save('data/inchikey_test.npy', inchikey_test)
