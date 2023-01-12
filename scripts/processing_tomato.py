# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 15:27:35 2022

@author: DELL
"""

import os
import numpy as np
import pandas as pd
import pubchempy as pcp

from tqdm import tqdm
from matchms import Spectrum
from matchms.importing import load_from_mgf
from matchms.exporting import save_as_mgf, save_as_msp
from core.identification import spectrum_processing
from core.msdial import load_MS_DIAL_Alginment, remove_duplicate

'''
ident_pos = pd.read_excel("E:/Data/tomato_data/MH2102-009-DR-TMT-pos-LS-FINAL-MSMS.xlsx").loc[:,['Compound Name', 'm/z (Expected)', 'RT']]
ident_neg = pd.read_excel("E:/Data/tomato_data/MH2102-009-TMT-neg-LS.xlsx").loc[:,['Compound Name', 'm/z (Expected)', 'RT']]

ident_pos['Ionization'] = 'POS'
ident_neg['Ionization'] = 'NEG'

def remove_uncertain(data):
    keep = []
    for i in data.index:
        mz = data.loc[i, 'm/z (Expected)']
        rt = data.loc[i, 'RT']
        if np.isnan(mz) or np.isnan(rt):
            continue
        
        check = np.logical_and(np.abs(data.loc[:, 'm/z (Expected)'].values - mz) < 0.05, 
                               np.abs(data.loc[:, 'RT'].values - rt) < 0.5)
        if len(np.where(check)[0]) > 1:
            continue
        else:
            keep.append(i)
            
    final = data.loc[keep,:]
    final = final.reset_index(drop = True)
    return final

ident_pos = remove_uncertain(ident_pos)
ident_neg = remove_uncertain(ident_neg)
tomato_identified = ident_pos.append(ident_neg, ignore_index = True)


exact_mass, inchikey, molecular_formula = [], [], []
for i in tqdm(tomato_identified.index):
    name = tomato_identified.loc[i, 'Compound Name']
    if name[-1] == ' ':
        name = name[:-1]
    compounds = pcp.get_compounds(name, 'name')

    m = np.nan
    k = np.nan
    f = np.nan
    if len(compounds) > 0:
        cmd = compounds[0]
        m = float(cmd.exact_mass)
        k = cmd.inchikey
        f = cmd.molecular_formula
    else:
        name = name.replace('-', '')
        compounds = pcp.get_compounds(name, 'name')
        if len(compounds) > 0:
            cmd = compounds[0]
            m = float(cmd.exact_mass)
            k = cmd.inchikey
            f = cmd.molecular_formula
    
    exact_mass.append(m)
    inchikey.append(k)
    molecular_formula.append(f)

tomato_identified['exact_mass'] = exact_mass
tomato_identified['inchikey'] = inchikey
tomato_identified['molecular_formula'] = molecular_formula
tomato_identified.to_csv('example/Tomato/tomato_identified.csv', index = False)
'''

pos_path = 'example/Tomato/tomato_positive_msdial.csv'
neg_path = 'example/Tomato/tomato_negative_msdial.csv'
pos_data = pd.read_csv(pos_path)
neg_data = pd.read_csv(neg_path)
pos_cols = list(pos_data.columns[32:])
neg_cols = list(neg_data.columns[32:])

spectrums_positive = load_MS_DIAL_Alginment('example/Tomato/tomato_positive_msdial.csv', sample_cols = pos_cols)
spectrums_negative = load_MS_DIAL_Alginment('example/Tomato/tomato_negative_msdial.csv', sample_cols = neg_cols)

spectrums_positive = [s for s in spectrums_positive if len(s.intensities[s.intensities > 0.05]) >= 3]
spectrums_positive = [s for s in spectrums_positive if s.get('parent_mass') is not None]
spectrums_negative = [s for s in spectrums_negative if len(s.intensities[s.intensities > 0.05]) >= 3]
spectrums_negative = [s for s in spectrums_negative if s.get('parent_mass') is not None]
spectrums_positive = [s.set('compound_name', 'Compound_{}'.format(i)) for i, s in enumerate(spectrums_positive)]
spectrums_negative = [s.set('compound_name', 'Compound_{}'.format(i)) for i, s in enumerate(spectrums_negative)]
save_as_mgf(spectrums_positive, 'example/Tomato/ms_ms_tomato_all_positive.mgf')
save_as_mgf(spectrums_negative, 'example/Tomato/ms_ms_tomato_all_negative.mgf')

spectrums = spectrums_positive + spectrums_negative
spectrums = remove_duplicate(spectrums)

identified_spectrums = []
tomato_identified = pd.read_csv('example/Tomato/tomato_identified.csv')
tomato_identified = tomato_identified.dropna()

for s in tqdm(spectrums):
    mz_diff = np.abs(tomato_identified['m/z (Expected)'].values - s.metadata['precursor_mz'])
    rt_diff = np.abs(tomato_identified['RT'].values - s.metadata['retention_time'])
    type_check = tomato_identified['Ionization'].values == s.metadata['ionmode'][:3].upper()
    wh = np.where(np.logical_and(mz_diff <= 0.01, rt_diff <= 0.3, type_check))[0]
    if len(wh) > 0:
        wh = wh[0]
        name = tomato_identified.loc[wh, 'Compound Name']
        inchikey = tomato_identified.loc[wh, 'inchikey']
        
        s = s.set('inchikey', inchikey)
        s = s.set('true_annotation', name)
        identified_spectrums.append(s)
identified_spectrums = [s.set('compound_name', 'Compound_{}'.format(i)) for i, s in enumerate(identified_spectrums)]


save_as_mgf(identified_spectrums, 'example/Tomato/ms_ms_tomato_identified.mgf')

from core.msdial import save_as_sirius
export_path = 'example/Tomato/msp'
save_as_sirius(identified_spectrums, export_path)

from core.msdial import save_as_msfinder
export_path = 'example/Tomato/msfinder'
save_as_msfinder(identified_spectrums, export_path)