# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 08:57:48 2022

@author: DELL
"""

import os
import pickle
from tqdm import tqdm
from matchms.filtering import add_losses


# positive
path_data = 'D:/All_MSDatabase'

outfile = os.path.join(path_data, 'GNPS_all/ALL_GNPS_220601_positive_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums = pickle.load(file)
spectrums = [add_losses(s) for s in tqdm(spectrums)]
pickle.dump(spectrums, open(outfile, "wb"))


outfile = os.path.join(path_data, 'In_House/ALL_Inhouse_positive_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums = pickle.load(file)
spectrums = [add_losses(s) for s in tqdm(spectrums)]
pickle.dump(spectrums, open(outfile, "wb"))


outfile = os.path.join(path_data, 'NIST2020/ALL_NIST20_positive_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums = pickle.load(file)
spectrums = [add_losses(s) for s in tqdm(spectrums)]
pickle.dump(spectrums, open(outfile, "wb"))


# negative
outfile = os.path.join(path_data, 'GNPS_all/ALL_GNPS_220601_negative_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums = pickle.load(file)
spectrums = [add_losses(s) for s in tqdm(spectrums)]
pickle.dump(spectrums, open(outfile, "wb"))


outfile = os.path.join(path_data, 'In_House/ALL_Inhouse_negative_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums = pickle.load(file)
spectrums = [add_losses(s) for s in tqdm(spectrums)]
pickle.dump(spectrums, open(outfile, "wb"))


outfile = os.path.join(path_data, 'NIST2020/ALL_NIST20_negative_cleaned.pickle')
with open(outfile, 'rb') as file:
    spectrums = pickle.load(file)
spectrums = [add_losses(s) for s in tqdm(spectrums)]
pickle.dump(spectrums, open(outfile, "wb"))


