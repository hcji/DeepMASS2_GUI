# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 09:16:02 2022

@author: DELL
"""


import os
import requests
import numpy as np
import pandas as pd

from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from matchms import Spectrum
from matchms import filtering as msfilters
from matchms.importing import load_from_mgf
from matchms.exporting import save_as_mgf, save_as_msp


ADDUCT_MASSES = {
    # POSITIVE
    '[M+H]+': 1.007276,
    '[M-H2O+H]+': 1.007276 - 18.010565, '[M+H-H2O]+': 1.007276 - 18.010565,
    '[M-2H2O+H]+': 1.007276 - 2 * 18.010565, '[M+H-2H2O]+': 1.007276 - 2 * 18.010565,
    '[M+Na]+': 22.98922,
    '[M]+': 0,
    '[M+NH4]+': 18.03383,
    '[M+H-NH3]+': 1.007276 - 17.026549,
    # NEGATIVE
    '[M-H]-': -1.007276,
    '[M+Cl]-': 34.969402,
    '[M+FA-H]-': 44.9982, '[M-H+FA]-': 44.9982,
    '[M]-': 0,
    '[M-2H]-': -1.007276
}


def spectrum_processing(s):
    """This is how one would typically design a desired pre- and post-
    processing pipeline."""
    s = msfilters.default_filters(s)
    s = msfilters.add_parent_mass(s)
    s = msfilters.normalize_intensities(s)
    s = msfilters.select_by_mz(s, mz_from=0, mz_to=2000)
    s = msfilters.add_losses(s, loss_mz_from=10.0, loss_mz_to=200.0)
    return s


def get_true_precursor_mz_from_mass(mass, precursor_type):
    """
    Calculate precursor mz based on (exact or monoisotopic) mass and precursor type
    :param mass: scalar, mass of a compound, e.g. monoisotopic or exact mass.
    :param precursor_type: string, precursor type, e.g. '[M+H]+'
    :return: scalar, ion mass / precursor mz
    """
    try:
        return mass + ADDUCT_MASSES[precursor_type]
    except KeyError:
        raise KeyError("Unsupported precursor-type '%s'." % precursor_type)


def get_adduct_from_mass_precursor(exactmass, precursor_mz, ion_mode, tolerence):
    if ion_mode == 'positive':
        k = [k for k in ADDUCT_MASSES.keys() if k[-1] == '+']
    else:
        k = [k for k in ADDUCT_MASSES.keys() if k[-1] == '-']
    adduct = ''
    diff = precursor_mz - exactmass
    for kk in k:
        if abs(diff - ADDUCT_MASSES[kk]) <= tolerence:
            adduct = kk
            break
    return adduct


challenge_ms = []

# CASMI 2022
files = os.listdir('d:/All_CASMI/data/CASMI2022_Priority')
challenge = pd.read_csv('d:/All_CASMI/solutions/casmi_2022_challenge_priority.csv')

spectrums = []
for f in tqdm(files):
    f_ = 'd:/All_CASMI/data/CASMI2022_Priority/{}'.format(f)
    mode = f.split('_')[2][:3]
    challenge_ = challenge[challenge['File'] == f.split('.')[0]]
    spectrums_ = [s for s in load_from_mgf(f_)]
    spectrums_ = [s.set('file', f.split('.')[0]) for s in spectrums_]
    spectrums += spectrums_


for i in tqdm(challenge.index):
    f_ = challenge.loc[i, 'File']
    rt_ = challenge.loc[i, 'RT [min]'] * 60
    mz_ = challenge.loc[i, 'Precursor m/z (Da)']
    intensity_ = 0
    mode = f_.split('_')[2][:3]
    
    name = 'challenge_{}'.format(len(challenge_ms))
    index = 'CASMI2022_Priority_{}'.format(i)
    smi = challenge.loc[i, 'SMILES']
    inchikey = challenge.loc[i, 'InChIKey']
    precursor_type = challenge.loc[i, 'Adduct'].replace(' ', '')
    mol = Chem.MolFromSmiles(smi)
    exactmass = AllChem.CalcExactMolWt(mol)
    formula = AllChem.CalcMolFormula(mol)
    precursor_mz = get_true_precursor_mz_from_mass(exactmass, precursor_type)
    
    s_new = None
    for s in spectrums:
        if s.metadata['file'] != f_:
            continue
        elif abs(s.metadata['precursor_mz'] - mz_) >= 0.01:
            continue
        elif abs(s.metadata['retention_time'] - rt_) >= 5:
            continue
        elif s.metadata['precursor_intensity'] <= intensity_:
            continue
        else:
            intensity_ = s.metadata['precursor_intensity']
            s_new = Spectrum(mz = s.mz, intensities = s.intensities, 
                             metadata = {'name': name,
                                         'precursor_mz': precursor_mz})
            s_new.set('challenge', index)
            s_new.set('formula', formula)
            s_new.set('smiles', smi)
            s_new.set('inchikey', inchikey)
            s_new.set('adduct', precursor_type)
            if mode == 'pos':
                s_new.set('ionmode', 'positive')
            else:
                s_new.set('ionmode', 'negative')
            s_new = spectrum_processing(s_new)
            
    if s_new:
        challenge_ms.append(s_new)


# CASMI 2016
files = os.listdir('d:/All_CASMI/data/CASMI2016_Test')
challenge = pd.read_csv('d:/All_CASMI/solutions/casmi_2016_challenge_test.csv')

for i, f in enumerate(tqdm(files)):
    f_ = 'd:/All_CASMI/data/CASMI2016_Test/{}'.format(f)
    spectrum = [s for s in load_from_mgf(f_)][0]
    
    name = 'challenge_{}'.format(len(challenge_ms))
    index = 'CASMI2016_{}'.format(challenge.loc[i, 'ChallengeName'])
    smi = challenge.loc[i, 'SMILES']
    inchikey = challenge.loc[i, 'INCHIKEY']
    mol = Chem.MolFromSmiles(smi)
    exactmass = AllChem.CalcExactMolWt(mol)
    formula = AllChem.CalcMolFormula(mol)
    precursor_mz = challenge.loc[i, 'PRECURSOR_MZ']
    ion_mode = challenge.loc[i, 'ION_MODE'].lower().replace(' ', '')
    
    adduct = get_adduct_from_mass_precursor(exactmass, precursor_mz, ion_mode, tolerence = 0.1)
    precursor_mz = get_true_precursor_mz_from_mass(exactmass, adduct)
    
    s_new = Spectrum(mz = spectrum.mz, intensities = spectrum.intensities, 
                     metadata = {'name': name,
                                 'precursor_mz': precursor_mz})
    s_new.set('challenge', index)
    s_new.set('formula', formula)
    s_new.set('smiles', smi)
    s_new.set('inchikey', inchikey)
    s_new.set('adduct', adduct)
    s_new.set('ionmode', ion_mode)
    s_new = spectrum_processing(s_new)
    challenge_ms.append(s_new)
    

# CASMI 2014
challenge = pd.read_csv('d:/All_CASMI/solutions/casmi_2014_challenge.csv')

for i in tqdm(challenge.index):
    link = 'http://casmi-contest.org/2014/Challenge2014/Challenge{}/{}_MSMS.txt'.format(i+1, i+1)
    assay = 'http://casmi-contest.org/2014/Challenge2014/Challenge{}/{}_Experimental_details.txt'.format(i+1, i+1)
    precursor = 'http://casmi-contest.org/2014/Challenge2014/Challenge{}/{}_MS.txt'.format(i+1, i+1)
    
    try:
        ms = requests.get(link).text.split('\n')
        ms = [s.replace('\r', '') for s in ms if s != '']
        ms = [s.replace(' ', '') for s in ms if s != '']
        mz = np.array([float(s.split('\t')[0]) for s in ms if s != ''])
    except:
        continue
    intensities = np.array([float(s.split('\t')[1]) for s in ms if s != ''])
    
    details = requests.get(assay).text
    if 'positive' in details.lower().split(' '):
        ion_mode = 'positive'
    else:
        ion_mode = 'negative'

    for r in requests.get(precursor).text.split('\n'):
        try:
            precursor_mz = float(r.split('\t')[0])
            continue
        except:
            pass
    
    name = 'challenge_{}'.format(len(challenge_ms))
    index = 'CASMI2016_{}'.format(challenge.loc[i, 'Challenge'])
    inchi = challenge.loc[i, 'InChI']
    mol = Chem.MolFromInchi(inchi)
    smi = Chem.MolToSmiles(mol)
    inchikey = challenge.loc[i, 'InChIkey']
    exactmass = AllChem.CalcExactMolWt(mol)
    formula = AllChem.CalcMolFormula(mol)
    
    adduct = get_adduct_from_mass_precursor(exactmass, precursor_mz, ion_mode, tolerence = 0.1)
    if adduct == '':
        if ion_mode == 'positive':
            adduct = '[M+H]+'
        else:
            adduct = '[M-H]-'
    precursor_mz = get_true_precursor_mz_from_mass(exactmass, adduct)
    
    s_new = Spectrum(mz = mz, intensities = intensities, 
                     metadata = {'name': name,
                                 'precursor_mz': precursor_mz})
    s_new.set('challenge', index)
    s_new.set('formula', formula)
    s_new.set('smiles', smi)
    s_new.set('inchikey', inchikey)
    s_new.set('adduct', adduct)
    s_new.set('ionmode', ion_mode)
    s_new = spectrum_processing(s_new)
    challenge_ms.append(s_new)


# save as mgf
save_as_mgf(list(challenge_ms), 'example/CASMI/all_casmi.mgf')

# save as msp individually
for i, s in enumerate(challenge_ms):
    path = 'example/CASMI/msp/challenge_{}.msp'.format(i)
    save_as_msp([challenge_ms[i]], path)
    
    with open(path) as msp:
        lines = msp.readlines()
        lines = [l.replace('_', '') for l in lines]
        lines = [l.replace('ADDUCT', 'PRECURSORTYPE') for l in lines]
    with open(path, 'w') as msp:
        msp.writelines(lines)

    # only for ms-finder
    path_msfinder = path.replace('/msp/', '/msfinder/')
        
    # exclude large compound, as ms-finder is very slow for them
    if s.metadata['parent_mass'] >= 850:
        continue
    with open(path_msfinder, 'w') as msp:
        msp.writelines(lines)
 
