# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 14:19:02 2022

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
from core.msdial import load_MS_DIAL_Alginment

def remove_duplicate(spectrums):
    new_spectrums = []
    rt, mz, iontype, intensities = [], [], [], []
    for s in tqdm(spectrums):
        [rt_, mz_, iontype_, intensity_, adduct_] = [s.metadata[k] for k in ['retention_time', 'precursor_mz', 'ionmode', 'precursor_intensity', 'adduct']]
        if adduct_ not in ['[M+H]+', '[M-H]-']:
            continue
        wh = np.logical_and( np.abs(np.array(rt) - rt_) < 18,
                             np.abs(np.array(mz) - mz_) < 0.01,
                             np.array([i == iontype_ for i in iontype]))
        wh = np.where(wh)[0]
        if len(wh) > 0:
            w = wh[0]
            if intensity_ >= intensities[w]:
                new_spectrums[w] = s
                intensities[w] = intensity_
            else:
                continue
        else:
            rt.append(rt_)
            mz.append(mz_)
            iontype.append(iontype_)
            intensities.append(intensity_)
            new_spectrums.append(spectrum_processing(s))
    return new_spectrums


spectrums = []
spectrums += load_MS_DIAL_Alginment('example/Plasma/ms_dial_positive.csv')
spectrums += load_MS_DIAL_Alginment('example/Plasma/ms_dial_negative.csv')
spectrums = [s for s in spectrums if len(s.intensities[s.intensities > 0.05]) >= 3]
spectrums = [s for s in spectrums if s.get('parent_mass') is not None]
spectrums = remove_duplicate(spectrums)
spectrums = [s.set('compound_name', 'Compound_{}'.format(i)) for i, s in enumerate(spectrums)]


'''
for f in tqdm(os.listdir("E:/Data/MTBLS263/positive/mgf")):
    f_ = 'E:/Data/MTBLS263/positive/mgf/{}'.format(f)
    mode = f.split('_')[2][:3]
    spectrums_ = [s for s in load_from_mgf(f_) if s.get('precursor_intensity') > 50000]
    spectrums_ = [s.set('ionmode', 'positive') for s in spectrums_]
    spectrums += spectrums_

    
for f in tqdm(os.listdir("E:/Data/MTBLS263/negative/mgf")):
    f_ = 'E:/Data/MTBLS263/negative/mgf/{}'.format(f)
    mode = f.split('_')[2][:3]
    spectrums_ = [s for s in load_from_mgf(f_) if s.get('precursor_intensity') > 50000]
    spectrums_ = [s.set('ionmode', 'negative') for s in spectrums_]
    spectrums += spectrums_


def preprocess_spectrums(spectrums):
    new_spectrums = []
    for i, s in enumerate(tqdm(spectrums)):
        name = 'unknown_{}'.format(i)
        if s.get('ionmode') == 'positive':
            adduct = '[M+H]+'
        else:
            adduct = '[M-H]-'
        
        s_new = Spectrum(mz = s.mz, intensities = s.intensities, 
                         metadata = {'name': name,
                                     'precursor_mz': s.get('precursor_mz')})
        s_new.set('adduct', adduct)
        s_new.set('retention_time', s.get('retention_time') / 60)
        s_new.set('ionmode', s.get('ionmode'))
        s_new.set('precursor_intensity', s.get('precursor_intensity'))
        s_new = spectrum_processing(s_new)
        new_spectrums.append(s_new)
    return new_spectrums


spectrums = remove_duplicate(spectrums)
spectrums = preprocess_spectrums(spectrums)
'''

'''
# retrieve information with name
pnas_identified = pd.read_csv('example/Plasma/pnas_identified.csv')
exact_mass, inchikey, molecular_formula = [], [], []
for i in tqdm(pnas_identified.index):
    name = pnas_identified.loc[i, 'Name']
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

pnas_identified['exact_mass'] = exact_mass
pnas_identified['inchikey'] = inchikey
pnas_identified['molecular_formula'] = molecular_formula
pnas_identified.to_csv('example/Plasma/pnas_identified.csv', index = False)
'''

# Manual correction
pnas_identified = pd.read_csv('example/Plasma/pnas_identified.csv')
pnas_identified['Ms Diff'] = pnas_identified['Detected m/z (A)'] - pnas_identified['exact_mass']

k = 0
for s in spectrums:
    mz_diff = np.abs(pnas_identified['Detected m/z (A)'].values - s.metadata['precursor_mz'])
    rt_diff = np.abs(pnas_identified['RT (min)'].values - s.metadata['retention_time'])
    type_check = pnas_identified['Ionization'].values == s.metadata['ionmode'][:3].upper()
    wh = np.where(np.logical_and(mz_diff <= 0.01, rt_diff <= 0.3, type_check))[0]
    if len(wh) > 0:
        wh = wh[0]
        name = pnas_identified.loc[wh, 'Name']
        inchikey = pnas_identified.loc[wh, 'inchikey']
        
        s.set('inchikey', inchikey)
        s.set('true_annotation', name)
        k += 1
save_as_mgf(spectrums, 'example/Plasma/ms_ms_plasma.mgf')


from core.msdial import save_as_sirius
export_path = 'example/Plasma/msp'
save_as_sirius(spectrums, export_path)

from core.msdial import save_as_msfinder
export_path = 'example/Plasma/msfinder'
save_as_msfinder(spectrums, export_path)