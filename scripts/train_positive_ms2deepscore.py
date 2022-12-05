# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 14:28:08 2022

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


tanimoto_scores_file = os.path.join(path_data, 'ALL_positive_tanimoto_scores.pickle')
with open(tanimoto_scores_file, 'rb') as file:
    tanimoto_scores_df = pickle.load(file)
    
tanimoto_scores_df = tanimoto_scores_df.fillna(0)


from ms2deepscore import SpectrumBinner
spectrum_binner = SpectrumBinner(1000, mz_min=10.0, mz_max=1000.0, peak_scaling=0.5, allowed_missing_percentage=99)
binned_spectrums = spectrum_binner.fit_transform(spectrums)


from ms2deepscore.data_generators import DataGeneratorAllSpectrums
dimension = len(spectrum_binner.known_bins)
data_generator = DataGeneratorAllSpectrums(binned_spectrums, tanimoto_scores_df, dim=dimension)


from tensorflow import keras
from ms2deepscore.models import SiameseModel


epochs = 30
model = SiameseModel(spectrum_binner, base_dims=(300, 300, 300), embedding_dim=200,
                     dropout_rate=0.2)
model.compile(loss='mse', optimizer=keras.optimizers.Adam(lr=0.001))
model.fit(data_generator, epochs = epochs)
model.save(os.path.join(path_data, 'MS2DeepScore_allGNPSpositive.hdf5'))

