# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 09:13:11 2022

@author: DELL
"""

import os
import numpy as np
import pandas as pd

from tqdm import tqdm
from matchms.importing import load_from_mgf
from matchms.exporting import save_as_mgf
import matchms.filtering as msfilters

challenge = pd.read_csv('d:/DeepMASS2/data/casmi_2022_challenge_priority.csv')
challenge_ms = np.repeat(None, len(challenge))
challenge_mode = np.repeat(None, len(challenge))
pred_formula = np.repeat(None, len(challenge))
files = os.listdir('d:/DeepMASS2/data/CASMI2022')

prdicted = pd.read_csv("E:\\Data\\Metabolomics Society 2022 - Workshop - UC Davis\\final_priority_Nothias.csv")
for i in challenge.index:
    cid = challenge.loc[i, 'Compound Number']
    w = np.where(prdicted['Compound Number'] == cid)[0]
    if len(w) == 0:
        continue
    else:
        fm = prdicted.loc[w[0], 'Molecular Formula']
        pred_formula[i] = fm
challenge['Predicted Formula'] = pred_formula
    

def spectrum_processing(s):
    """This is how one would typically design a desired pre- and post-
    processing pipeline."""
    s = msfilters.default_filters(s)
    s = msfilters.add_parent_mass(s)
    s = msfilters.normalize_intensities(s)
    s = msfilters.select_by_mz(s, mz_from=0, mz_to=1000)
    s = msfilters.add_losses(s, loss_mz_from=10.0, loss_mz_to=200.0)
    return s

for f in tqdm(files):
    f_ = 'd:/DeepMASS2/data/CASMI2022/{}'.format(f)
    mode = f.split('_')[2][:3]
    challenge_ = challenge[challenge['File'] == f.split('.')[0]]
    spectrums = load_from_mgf(f_)
    for s in spectrums:
        if np.min(np.abs(s.metadata['precursor_mz'] - challenge_['Precursor m/z (Da)'])) < 0.01:
            w = np.argmin(np.abs(s.metadata['precursor_mz'] - challenge_['Precursor m/z (Da)']))
            if (challenge_['RT [min]'].values[w] - s.metadata['retention_time']) < 0.5:
                i = challenge_.index[w]
                formula = challenge_.loc[i, 'Predicted Formula']
                precursor_type = challenge_.loc[i, 'Adduct']
                s.set('formula', formula)
                s.set('precursor_type', precursor_type)
                if mode == 'pos':
                    s.set('ionmode', 'positive')
                else:
                    s.set('ionmode', 'negative')
                if challenge_ms[i] is None:
                    challenge_ms[i] = spectrum_processing(s)
                    challenge_mode[i] = mode
                else:
                    if 'precursor_intensity' not in list(s.metadata.keys()):
                        continue
                    if challenge_ms[i].metadata['precursor_intensity'] < s.metadata['precursor_intensity']:
                        challenge_ms[i] = spectrum_processing(s)
                        challenge_mode[i] = mode
                    

save_as_mgf(list(challenge_ms), 'example/CASMI_2022.mgf')
